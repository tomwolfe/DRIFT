"""Solver fallbacks and seed propagation.

Two properties that are easy to believe and easy to lose:

  * the SOLVER FALLBACK: when no COBRA-compatible LP solver is available the
    engine must degrade to the headless proxy (and refuse to run in strict
    mode) rather than pretend it solved an LP;
  * the SEED PROPAGATION: a seeded run must be reproducible, and the seed
    must reach BOTH the integrator's own generator and the Numba global
    generator the jitted drift uses.

The fallback is asserted by MEASURING the solver registry, not by reading
the code: the test asks optlang what it has, and then checks that the
engine's notion of "no solver" agrees with it.  Where the registry is
emptied, it is emptied by patching what optlang REPORTS, so the engine
reads exactly the measurement the test asserts on -- no hand-set
``available_solvers``/``headless`` attributes, no monkeypatched
``_check_solver``.

Where an environment genuinely cannot run a case (no LP solver installed),
the test SKIPS with the measured registry in the reason -- it never passes
vacuously.

Note on the two cases that were previously asserted more weakly: the
"missing solver" case used to hand-set ``available_solvers``/``headless``
on a solver that still had a model bound (so ``optimize()`` returned
"optimal" and the proxy was never exercised), and the strict case used to
monkeypatch ``_check_solver`` away.  Both are now driven by the measured
registry, which is strictly stronger: the engine really reads an empty
registry and really refuses.
"""

from __future__ import annotations

from typing import List

import numpy as np
import pytest
from numba import njit

from drift.binding import BindingEngine
from drift.engine import SimulationEngine
from drift.metabolic import DFBASolver, MetabolicBridge
from drift.signaling import StochasticIntegrator, seed_numba
from drift.topology import get_default_topology

try:  # cobra drags optlang in, but the registry is what we measure, so be explicit
    import optlang
except ImportError:  # pragma: no cover - optlang missing => registry is empty
    optlang = None  # type: ignore[assignment]


@njit
def _numba_draw() -> float:
    """Draw from NUMBA's generator.

    Numba's RNG is a separate state from NumPy's, which is exactly why the
    seeding call lives in a jitted function: the jitted drift draws its
    noise from Numba's generator.  Tests therefore draw from a jitted
    function, never from ``np.random`` (a different generator).
    """
    return np.random.rand()


def _installed_solvers() -> List[str]:
    """MEASURE the solver registry: the solvers optlang actually reports."""
    if optlang is None:
        return []
    return [name for name, available in optlang.available_solvers.items() if available]


def _empty_registry(monkeypatch) -> None:
    """Make the registry MEASURE as empty, so the engine must degrade."""
    if optlang is None:
        # Without optlang the engine's own check already reports no solver.
        return
    monkeypatch.setattr(optlang, "available_solvers", {}, raising=False)


def _make_engine(noise_scale: float = 0.01) -> SimulationEngine:
    """A stock engine: textbook metabolic model, default signaling topology."""
    integrator = StochasticIntegrator(
        dt=0.05, noise_scale=noise_scale, topology=get_default_topology())
    solver = DFBASolver(model_name="textbook")
    bridge = MetabolicBridge()
    binding = BindingEngine(targets=1.0, hill_coefficient=1.0)
    return SimulationEngine(integrator, solver, bridge, binding)


# ---------------------------------------------------------------------------
# SOLVER FALLBACK
# ---------------------------------------------------------------------------


def test_solver_registry_agrees_with_optlang() -> None:
    """The engine's solver list IS the measured registry, not a guess."""
    installed = _installed_solvers()
    solver = DFBASolver(model_name="textbook")
    assert set(solver.available_solvers) == set(installed)
    assert solver.headless == (not solver.available_solvers)


def test_missing_solver_degrades_to_the_headless_proxy(monkeypatch) -> None:
    """An empty measured registry must degrade to the qualitative proxy.

    The registry is emptied by patching what optlang reports, so
    ``_check_solver`` really reads an empty registry.  The result must be
    the proxy -- never an "optimal" verdict from a solver that is not there.
    """
    _empty_registry(monkeypatch)
    assert _installed_solvers() == []

    solver = DFBASolver(model_name="textbook")
    assert solver.headless is True
    assert solver.available_solvers == []
    assert solver.model is None, "a model was bound although no solver exists"

    result = solver.optimize({"R_EX_lactate_e": -1.0})
    assert result["status"] == "proxied", result
    assert "QUALITATIVE" in result["diagnostic"].upper(), result
    assert result["objective_value"] > 0.0, result


def test_the_engine_runs_headless_when_the_registry_is_empty(monkeypatch) -> None:
    """A full run with no solver is headless and proxied at every step."""
    _empty_registry(monkeypatch)

    engine = _make_engine()
    result = engine.run(steps=3, drug_concentration=0.5, seed=5)

    assert result["headless"] is True
    assert result["status"] == ["proxied", "proxied", "proxied"], result["status"]
    assert "optimal" not in result["status"]


def test_strict_mode_refuses_to_run_without_a_solver(monkeypatch) -> None:
    """Strict mode must REFUSE, not degrade, when no solver is available."""
    _empty_registry(monkeypatch)
    assert _installed_solvers() == []

    with pytest.raises(RuntimeError, match="strict"):
        DFBASolver(model_name="textbook", strict=True)

    # The contrast, on the very same measurement: without strict the engine
    # degrades to the proxy instead of refusing.
    solver = DFBASolver(model_name="textbook", strict=False)
    assert solver.headless is True
    assert solver.optimize({"R_EX_lactate_e": -1.0})["status"] == "proxied"


def test_strict_mode_refuses_at_solve_time(monkeypatch) -> None:
    """The second guard: strict mode must not quietly proxy when solving."""
    _empty_registry(monkeypatch)

    solver = DFBASolver(model_name="textbook", strict=False)
    assert solver.headless is True
    assert solver.model is None

    solver.strict = True
    with pytest.raises(RuntimeError, match="strict"):
        solver.optimize({"R_EX_lactate_e": -1.0})


def test_a_solver_makes_the_engine_non_headless() -> None:
    """With a solver available the engine is NOT headless.

    The GLPK case is the one the fallback exists for: optlang's registry
    lists GLPK when swiglpk is installed, and the engine must then use a
    real LP solver instead of the proxy.  When GLPK is not installed the
    test is skipped WITH A REASON rather than passing vacuously.
    """
    installed = [name.lower() for name in _installed_solvers()]
    if not installed:
        pytest.skip(f"no solver installed (measured: {installed or 'none'})")
    target = installed[0]
    solver = DFBASolver(model_name="textbook")
    assert target in [name.lower() for name in solver.available_solvers]
    assert solver.headless is False
    assert solver.model is not None, "a solver is available but no model was bound"
    interface = solver.model.solver.interface.__name__.lower()
    assert target in interface, (
        f"the model is bound to {interface}, not the measured {target} solver")


def test_strict_mode_runs_when_a_solver_is_available() -> None:
    """Strict mode is only a refusal when there is genuinely no solver."""
    installed = _installed_solvers()
    if not installed:
        pytest.skip(f"no COBRA-compatible LP solver is installed in this "
                    f"environment (measured registry: {installed or 'empty'})")

    solver = DFBASolver(model_name="textbook", strict=True)
    assert solver.headless is False
    assert solver.available_solvers

    result = solver.optimize({"R_EX_lactate_e": -1.0})
    assert result["status"] != "proxied", (
        "strict mode proxied although a real solver is available")
    assert result["status"] in ("optimal", "infeasible"), result


# ---------------------------------------------------------------------------
# SEED PROPAGATION
# ---------------------------------------------------------------------------


def test_numba_seed_is_propagated() -> None:
    """``seed_numba`` really seeds NUMBA's generator.

    Numba's RNG is a separate state from NumPy's, which is exactly why the
    seeding call lives in a jitted function: the jitted drift draws its
    noise from Numba's generator, so only Numba's state matters for the
    simulation.  The test therefore draws from a jitted function, not from
    ``np.random`` (which is a different generator).
    """
    seed_numba(1234)
    first = _numba_draw()
    seed_numba(1234)
    second = _numba_draw()
    assert first == second, (
        "seed_numba did not reseed Numba's generator: the jitted drift's "
        "noise is not reproducible")


def test_the_engine_seeds_numba() -> None:
    """``run(seed=...)`` must reseed NUMBA's generator.

    The jitted drift draws its noise from Numba's generator, so a run whose
    seed never reaches Numba is not reproducible in the only sense that
    matters.  Numba's state is deliberately perturbed before the run: a
    stale stream would otherwise produce the same number twice and let a
    dropped ``seed_numba`` call hide.
    """
    engine = _make_engine()

    seed_numba(4242)
    expected = _numba_draw()
    _numba_draw()  # perturb Numba's global state

    engine.run(steps=4, drug_concentration=0.5, seed=4242)
    assert _numba_draw() == expected, (
        "the seed did not reach Numba's generator: the jitted drift's noise "
        "is not seeded by run(seed=...)")

    # A different seed must reach Numba's generator as well.
    seed_numba(777)
    other = _numba_draw()
    _numba_draw()  # perturb again

    engine.run(steps=4, drug_concentration=0.5, seed=777)
    assert _numba_draw() == other
    assert other != expected, "two different seeds gave the same Numba stream"


def test_the_integrator_seeds_both_generators() -> None:
    """``set_seed`` must seed NumPy's generator AND Numba's."""
    integrator = StochasticIntegrator(
        dt=0.05, noise_scale=0.01, topology=get_default_topology())
    integrator.set_seed(42)

    # NumPy's side: the integrator's own generator is the seeded stream.
    assert integrator.rng.random() == np.random.default_rng(42).random()

    # Numba's side: the jitted drift's generator was reseeded too.
    seed_numba(42)
    expected = _numba_draw()
    _numba_draw()  # perturb Numba's global state

    integrator.set_seed(42)
    assert _numba_draw() == expected, (
        "set_seed did not reach Numba's generator")


def test_the_engine_seeds_the_engine_and_integrator_generators() -> None:
    """``run(seed=...)`` seeds the engine's generator and the integrator's."""
    engine = _make_engine()
    species = len(engine.integrator.topology.species)

    engine.run(steps=1, drug_concentration=0.5, seed=42, sub_steps=1)

    # The engine's generator is the seeded stream, advanced by exactly the
    # noise draw the integrator took during the run.
    expected = np.random.default_rng(42)
    expected.standard_normal(species)
    assert engine.rng.random() == expected.random(), (
        "the seed did not reach the engine's generator")
    assert engine.integrator.rng is engine.rng, (
        "the integrator is not using the engine's seeded generator")


def test_a_seeded_simulation_is_reproducible() -> None:
    """The same seed gives the same trajectory, through the engine."""
    def _run():
        engine = _make_engine()
        return engine.run(steps=8, drug_concentration=0.5, seed=2024)

    first = _run()
    second = _run()
    assert np.asarray(first["signaling"]).shape == np.asarray(second["signaling"]).shape
    assert np.array_equal(np.asarray(first["signaling"]),
                          np.asarray(second["signaling"])), (
        "a seeded simulation is not reproducible: the seed did not reach "
        "every generator the integrator uses")
    assert np.array_equal(np.asarray(first["growth"]),
                          np.asarray(second["growth"]))
    assert first["status"] == second["status"]


def test_a_different_seed_gives_a_different_trajectory() -> None:
    """A seed that is not propagated would make this fail."""
    def _run(seed: int):
        engine = _make_engine(noise_scale=0.05)
        return engine.run(steps=8, drug_concentration=0.5, seed=seed)

    a = _run(1)
    b = _run(2)
    assert not np.allclose(np.asarray(a["signaling"]),
                           np.asarray(b["signaling"])), (
        "two different seeds gave the same trajectory: the seed is not "
        "reaching the integrator's noise")
    assert not np.allclose(np.asarray(a["growth"]),
                           np.asarray(b["growth"]))
