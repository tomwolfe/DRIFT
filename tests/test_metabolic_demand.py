"""A metabolic DEMAND, and the cell death it really causes.

A capacity mapping caps what a reaction may do (the historical
behaviour, and it must not change).  A requirement mapping states what a
reaction MUST carry: the bound is positive, the flux has to reach it, and
a demand above what the network can produce makes the LP INFEASIBLE.

The cell-death flag must be reached that way -- the solver's own verdict,
observed by the engine -- and never by setting the flag.  These tests
pin both halves:

  * the capacity path is byte-identical to what it was;
  * a demand above the network's capacity produces a genuinely infeasible
    LP, and the ENGINE notices and sets the flag and zeroes the growth;
  * the same wiring with a feasible demand does NOT kill the cell, so the
    flag tracks the infeasibility and not the plumbing.

Nothing here monkeypatches the solver: the status the engine reacts to
is the status the LP solver returned.
"""

from __future__ import annotations

import numpy as np
import pytest

from drift.binding import BindingEngine
from drift.engine import SimulationEngine
from drift.metabolic import DFBASolver, MetabolicBridge, constraint_bound
from drift.signaling import StochasticIntegrator
from drift.topology import get_default_topology

BIOMASS_REACTION = "Biomass_Ecoli_core"
INFEASIBLE_DEMAND = 5.0
FEASIBLE_DEMAND = 0.1


def _bridge(demand: float | None = None) -> MetabolicBridge:
    """The default bridge, plus a biomass REQUIREMENT when asked for one."""
    mappings = [
        {
            "protein_name": "mTOR",
            "reaction_id": "EX_glc__D_e",
            "influence": "positive",
            "base_vmax": 10.0,
            "basal_scaling": 0.1,
            "max_scaling": 0.9,
        },
    ]
    if demand is not None:
        mappings.append({
            "protein_name": "mTOR",
            "reaction_id": BIOMASS_REACTION,
            "constraint_kind": "requirement",
            "influence": "positive",
            "base_vmax": demand,
            "basal_scaling": 1.0,
            "max_scaling": 0.0,
        })
    return MetabolicBridge(mappings=mappings,
                           species_names=["PI3K", "AKT", "mTOR"])


def _engine(bridge: MetabolicBridge) -> SimulationEngine:
    integrator = StochasticIntegrator(
        dt=0.05, noise_scale=0.01, topology=get_default_topology())
    solver = DFBASolver(model_name="textbook")
    binding = BindingEngine(targets=0.5, hill_coefficient=1.5)
    return SimulationEngine(integrator, solver, bridge, binding)


# ---------------------------------------------------------------------------
# The bridge
# ---------------------------------------------------------------------------


def test_the_default_mapping_is_still_a_capacity() -> None:
    """The historical behaviour: a capacity is a NEGATIVE lower bound."""
    bridge = _bridge()
    state = np.asarray(get_default_topology().get_initial_state())
    constraints = bridge.get_constraints(state)
    assert constraints["EX_glc__D_e"] < 0.0
    # The bound is exactly -(base_vmax * scaling), as it always was.
    scaling = bridge.get_constraints_with_scalings(state)[1]["EX_glc__D_e"]
    assert constraints["EX_glc__D_e"] == pytest.approx(-10.0 * scaling)


def test_a_requirement_states_a_demand() -> None:
    """A requirement is a POSITIVE lower bound: the flux must reach it."""
    bridge = _bridge(demand=INFEASIBLE_DEMAND)
    state = np.asarray(get_default_topology().get_initial_state())
    constraints, scalings = bridge.get_constraints_with_scalings(state)
    assert constraints[BIOMASS_REACTION] > 0.0
    assert constraints[BIOMASS_REACTION] == pytest.approx(INFEASIBLE_DEMAND)
    # The capacity mapping around it is untouched.
    assert constraints["EX_glc__D_e"] < 0.0
    assert set(scalings) == {"EX_glc__D_e", BIOMASS_REACTION}


def test_the_two_constraint_kinds_are_what_they_claim() -> None:
    """The helper is the single place the two kinds are decided."""
    assert constraint_bound({}, 10.0, 0.5) == pytest.approx(-5.0)
    assert constraint_bound({"constraint_kind": "capacity"}, 10.0, 0.5) == \
        pytest.approx(-5.0)
    assert constraint_bound({"constraint_kind": "requirement"}, 10.0, 0.5) == \
        pytest.approx(5.0)


def test_an_unknown_constraint_kind_is_refused() -> None:
    """A typo must not turn a demand into a cap."""
    with pytest.raises(ValueError):
        constraint_bound({"constraint_kind": "demmand"}, 1.0, 1.0)
    bridge = MetabolicBridge(
        mappings=[{
            "protein_name": "mTOR",
            "reaction_id": "EX_glc__D_e",
            "constraint_kind": "demmand",
            "base_vmax": 10.0,
        }],
        species_names=["PI3K", "AKT", "mTOR"],
    )
    state = np.asarray(get_default_topology().get_initial_state())
    with pytest.raises(ValueError):
        bridge.get_constraints(state)


# ---------------------------------------------------------------------------
# The cell death
# ---------------------------------------------------------------------------


def test_an_infeasible_demand_makes_the_engine_report_the_death() -> None:
    """The LP is infeasible, and the ENGINE notices and sets the flag."""
    bridge = _bridge(demand=INFEASIBLE_DEMAND)
    engine = _engine(bridge)

    # The problem itself, solved by the LP solver.
    state = engine.integrator.topology.get_initial_state()
    constraints, scalings = bridge.get_constraints_with_scalings(state)
    direct = engine.solver.optimize(constraints, scalings=scalings)
    assert direct["status"] == "infeasible"
    assert BIOMASS_REACTION in direct["bottlenecks"]

    # The engine's own run over the same problem.
    history = engine.run(steps=4, drug_concentration=1.0, seed=3)
    assert history["cell_death"] is True
    assert history["death_step"] == 0
    assert BIOMASS_REACTION in (history["death_cause"] or "")
    assert list(history["growth"]) == [0.0, 0.0, 0.0, 0.0]
    assert history["status"][0] == "infeasible"


def test_a_feasible_demand_does_not_kill_the_cell() -> None:
    """The control: the same wiring, a satisfiable demand, no death."""
    bridge = _bridge(demand=FEASIBLE_DEMAND)
    engine = _engine(bridge)
    state = engine.integrator.topology.get_initial_state()
    constraints, scalings = bridge.get_constraints_with_scalings(state)
    direct = engine.solver.optimize(constraints, scalings=scalings)
    assert direct["status"] == "optimal"

    history = engine.run(steps=4, drug_concentration=1.0, seed=3)
    assert history["cell_death"] is False
    assert history["death_step"] is None
    assert all(growth > 0.0 for growth in history["growth"])
