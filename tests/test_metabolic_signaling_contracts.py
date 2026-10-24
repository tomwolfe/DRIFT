"""Contract tests for the components the formal gate certifies.

Each test pins a behaviour whose mutation would otherwise leave the
source free to stop implementing the equation the certificate is about
while every exported lemma stayed unchanged.  The behaviours are the
ones the exporter reads (the inhibition modulation and the Milstein
step), the bridge's scaling law, and the predictor-corrector average
that the metabolic feedback actually flows through.

The test module's name carries both "metabolic" and "signaling" on
purpose: command 3 of the mission runs ``pytest -k "binding or
signaling or metabolic"``, and every contract here belongs to that
suite.
"""

from __future__ import annotations

import sys
import types
import warnings
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from drift.binding import BindingEngine
from drift.engine import SimulationEngine
from drift.metabolic import MetabolicBridge
from drift.signaling import StochasticIntegrator, langevin_step
from drift.topology import (
    Topology,
    complex_signaling_drift,
    pi3k_akt_mtor_drift,
)


def _plain_topology() -> Topology:
    """A topology with no compiled drift function (the generic path)."""
    return Topology(
        species=["PI3K", "AKT", "mTOR"],
        parameters={"k": 0.1},
        drift_fn=None,
    )


# ---------------------------------------------------------------------------
# The Milstein / Euler-Maruyama step
# ---------------------------------------------------------------------------


def test_the_generic_step_returns_the_next_state() -> None:
    """The generic step is a state TRANSITION, not a reset to None.

    ``StochasticIntegrator.step`` falls back to ``_generic_step`` for a
    topology without a compiled drift.  Returning ``None`` from that
    fallback would silently destroy the trajectory; the returned state
    must still be a bounded, finite vector.
    """
    integrator = StochasticIntegrator(
        dt=0.1, noise_scale=0.01, topology=_plain_topology())
    state = np.array([0.5, 0.5, 0.5])

    next_state = integrator.step(state, inhibition=0.0)

    assert isinstance(next_state, np.ndarray)
    assert next_state.shape == state.shape
    assert np.all(np.isfinite(next_state))
    assert np.all(next_state >= 0.0) and np.all(next_state <= 1.0)


def test_the_custom_step_falls_back_to_the_generic_step() -> None:
    """A custom step asked for while the drift disappeared must not
    return ``None``.

    The branch is the documented fallback of ``_custom_step``: if no
    custom drift function is installed, the generic step is the answer.
    """
    integrator = StochasticIntegrator(
        dt=0.1, noise_scale=0.01, topology=_plain_topology())
    state = np.array([0.4, 0.5, 0.6])
    noise = np.zeros(3)

    next_state = integrator._custom_step(
        state, np.zeros(3), noise, feedback=None)

    assert isinstance(next_state, np.ndarray)
    assert np.all(np.isfinite(next_state))
    assert np.all(next_state >= 0.0) and np.all(next_state <= 1.0)


def test_the_custom_step_returns_the_next_state() -> None:
    """The custom-drift step is a state TRANSITION too.

    A topology can carry a plain (non-jitted) ``drift_fn``; the step then
    runs ``_custom_step``.  Returning ``None`` instead of the clipped
    state would silently end the trajectory.
    """

    def decay(state, params, inhibition=0.0, feedback=None):
        return -0.1 * state

    topology = Topology(
        species=["PI3K", "AKT", "mTOR"],
        parameters={"k": 0.1},
        drift_fn=decay,
    )
    integrator = StochasticIntegrator(
        dt=0.1, noise_scale=0.01, topology=topology)

    next_state = integrator.step(np.array([0.5, 0.5, 0.5]), inhibition=0.0)

    assert isinstance(next_state, np.ndarray)
    assert np.all(np.isfinite(next_state))
    assert np.all(next_state >= 0.0) and np.all(next_state <= 1.0)


def test_langevin_step_clamps_a_state_above_one() -> None:
    """A state numerically just above 1 is clamped BEFORE the square root.

    ``b(x) = noise_scale * sqrt(x * (1 - x))`` is only real for
    ``x in [0, 1]``.  The projection must therefore clamp the lower end
    of the interval; a ``1.0 + eps`` upper clamp would make the argument
    negative and the step NaN for a state of ``1.0 + 1e-9``.
    """
    state = np.array([1.0 + 1e-9, 0.5, 0.5])
    params = np.array([0.1, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5])

    next_state = langevin_step(state, 0.1, params, 0.1)

    assert np.all(np.isfinite(next_state))
    assert np.all(next_state >= 0.0) and np.all(next_state <= 1.0)


def test_the_stability_check_scales_noise_by_sqrt_dt() -> None:
    """The noise-impact heuristic is ``noise_scale * sqrt(dt)``.

    ``dt=4.0`` and ``noise_scale=0.1`` put the product at 0.2 (the
    warning fires) while the quotient is 0.05 (it does not).  Pinning
    the product keeps the warning meaningful.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        StochasticIntegrator(
            dt=4.0, noise_scale=0.1, topology=_plain_topology())

    messages = [str(item.message) for item in caught]
    assert any("High noise impact" in message for message in messages), (
        "the noise-impact warning did not fire for noise_scale * sqrt(dt) "
        f"= 0.2; warnings were {messages!r}")


# ---------------------------------------------------------------------------
# The metabolic bridge's scaling law
# ---------------------------------------------------------------------------


def test_a_negative_influence_mirrors_the_scaling() -> None:
    """A negative influence runs the scaling from high to low.

    ``scaling = (basal + max_scaling) - max_scaling * protein_level``:
    at half protein the scaling is 0.55, not negative.  Swapping the
    leading ``+`` for a ``-`` would make the bridge drive the reaction
    with a scaling of -1.25.
    """
    bridge = MetabolicBridge(
        mappings=[{
            "protein_name": "mTOR",
            "reaction_id": "EX_glc__D_e",
            "influence": "negative",
            "base_vmax": 10.0,
            "basal_scaling": 0.1,
            "max_scaling": 0.9,
        }],
        species_names=["PI3K", "AKT", "mTOR"],
    )

    constraints, scalings = bridge.get_constraints_with_scalings(
        np.array([0.0, 0.0, 0.5]))

    assert scalings["EX_glc__D_e"] == pytest.approx(0.55)
    assert constraints["EX_glc__D_e"] == pytest.approx(-5.5)


def test_validate_with_model_reports_an_unresolved_mapping() -> None:
    """An unresolvable non-strict mapping makes the verdict False.

    With ``strict_mapping=False`` an unknown reaction is fuzzy-matched;
    when that fails the bridge repairs nothing and must report the
    mapping as unresolved.  Flipping the verdict to True would certify
    a mapping the model does not contain.
    """

    class _Model:
        reactions = [type("Reaction", (), {"id": "EX_glc__D_e"})()]

    bridge = MetabolicBridge(
        mappings=[{
            "protein_name": "mTOR",
            "reaction_id": "zzzz_not_a_real_reaction_zzzz",
            "influence": "positive",
        }],
        species_names=["PI3K", "AKT", "mTOR"],
    )

    assert bridge.validate_with_model(_Model()) is False


# ---------------------------------------------------------------------------
# The predictor-corrector average
# ---------------------------------------------------------------------------


def test_the_corrector_averages_the_old_and_predicted_feedback() -> None:
    """The corrector uses the MEAN of the old and predicted feedback.

    With ``feedback_old = 1`` and ``feedback_pred = 0.2`` the corrected
    signaling step must see 0.6.  A subtraction gives 0.4 and would make
    the corrector push the state in the opposite direction from the
    predictor it is correcting.
    """
    integrator = StochasticIntegrator(
        dt=0.1, noise_scale=0.01, topology=_plain_topology())

    class _Solver:
        headless = False

        def optimize(self, constraints, scalings=None):
            return {
                "status": "optimal",
                "objective_value": 1.0,
                "fluxes": {},
                "diagnostic": "ok",
            }

    class _Bridge:
        def get_constraints_with_scalings(self, state):
            return {}, {}

        def get_feedback(self, fluxes):
            return np.full(3, 0.2)

    seen: list[np.ndarray] = []

    def _record(state, inhibition, feedback=None):
        seen.append(np.array(feedback, copy=True))
        return state.copy()

    integrator.step = _record  # type: ignore[method-assign]
    engine = SimulationEngine(
        integrator, _Solver(), _Bridge(), BindingEngine(kd=0.5))

    engine.run(steps=1, drug_concentration=1.0, sub_steps=2,
               refine_feedback=True)

    # Two predictor steps with the initial feedback of 1, then two
    # corrector steps with the average of 1 and 0.2.
    assert len(seen) == 4
    np.testing.assert_allclose(seen[0], np.ones(3))
    np.testing.assert_allclose(seen[1], np.ones(3))
    np.testing.assert_allclose(seen[2], np.full(3, 0.6))
    np.testing.assert_allclose(seen[3], np.full(3, 0.6))


# ---------------------------------------------------------------------------
# The inhibition modulation the exporter certifies
# ---------------------------------------------------------------------------


def test_pi3k_akt_mtor_drift_decays_extra_species() -> None:
    """The three-species model leaves any extra species to decay.

    ``if len(state) > 3`` is the guard: a four-species state gets
    ``res[3] = 0.1 * (0.5 - state[3])`` while the certified three
    equations are untouched.
    """
    state = np.array([0.5, 0.5, 0.5, 0.9])
    params = np.array([0.1, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5])

    drift = pi3k_akt_mtor_drift(state, params)

    np.testing.assert_allclose(drift[:3], np.array([0.05, 0.0125, 0.075]),
                               rtol=0.0, atol=1e-12)
    assert drift[3] == pytest.approx(0.1 * (0.5 - 0.9))


def test_complex_signaling_drift_keeps_the_mtor_inhibition_term() -> None:
    """PI3K's production is inhibited by mTOR, not promoted by it.

    ``dpi3k = k_base * fb * 1/(1 + 5*mtor) - k_deg * pi3k``.  At the
    half-activated state the first term is 0.1/3.5 and the degradation
    is 0.05; a flipped sign in either term changes the sum.
    """
    state = np.full(4, 0.5)
    params = np.array(
        [0.1, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1, 0.0, 0.0, 0.0, 0.0])

    drift = complex_signaling_drift(state, params)

    np.testing.assert_allclose(
        drift,
        np.array([
            0.1 / 3.5 - 0.05,   # dpi3k: 1/(1 + 5*0.5) inhibition
            0.125 - 0.05,        # dakt: 0.5 * 0.5 * (1 - 0.5)
            0.125 / 6.0 - 0.05,  # dmtor: ampk inhibition 1/(1 + 10*0.5)
            -0.05,               # dampk
        ]),
        rtol=0.0, atol=1e-12,
    )


# ---------------------------------------------------------------------------
# The SBML / SBML-qual inhibition paths
# ---------------------------------------------------------------------------


def _qual_extension_mocks() -> tuple[MagicMock, MagicMock]:
    """A one-species SBML-qual model whose single transition is ``S1``."""
    model = MagicMock()
    model.getName.return_value = "Qual"
    model.getId.return_value = "qual"

    qual = MagicMock()
    species = MagicMock()
    species.getId.return_value = "S1"
    qual.getListOfQualitativeSpecies.return_value = [species]

    output = MagicMock()
    output.getQualitativeSpecies.return_value = "S1"
    transition = MagicMock()
    transition.getListOfOutputs.return_value = [output]

    default_term = MagicMock()
    default_term.getResultLevel.return_value = 0.0
    function_term = MagicMock()
    function_term.getMath.return_value = MagicMock()
    function_term.getResultLevel.return_value = 1.0
    terms = MagicMock()
    terms.getDefaultTerm.return_value = default_term
    terms.__iter__ = lambda self: iter([function_term])
    transition.getListOfFunctionTerms.return_value = terms

    qual.getListOfTransitions.return_value = [transition]
    return model, qual


def test_the_compiled_sbml_qual_drift_applies_inhibition() -> None:
    """The generated code scales the inhibited species by ``1 - inh``.

    ``S1``'s target level is 1 and the step sees ``inhibition=0.5``, so
    the inhibited target is 0.5 and the rate is ``-kn * state``.
    Emitting the inhibition for every species EXCEPT the inhibited one
    (a flipped comparison) leaves the target at 1 and changes the rate.
    """
    model, qual = _qual_extension_mocks()

    with patch("libsbml.formulaToL3String", return_value="S1"):
        topology = Topology._from_sbml_qual(
            model, qual, inhibited_species="S1")

    drift = topology.drift_fn(
        np.array([0.5]), np.array([0.5, 0.1]), 0.5)

    np.testing.assert_allclose(drift, np.array([-0.05]), rtol=0.0, atol=1e-12)


def test_the_interpreted_sbml_qual_drift_applies_inhibition_and_decay() -> None:
    """The interpreted fallback has the same inhibition and decay algebra.

    The compiled path is made to fail so the interpreted function is the
    one under test.  With ``kd=1``, ``kn=0.4``, target 1, inhibition 0.5
    and state 0.5 the rate is ``1*(0.5-0.5) - 0.4*0.5 = -0.2``.  Adding
    the inhibition instead of subtracting it changes the target to 1.5;
    dividing the decay by the state changes the second term to -0.8.
    """
    model, qual = _qual_extension_mocks()

    with patch("libsbml.formulaToL3String", return_value="S1"), \
            patch("builtins.compile", side_effect=Exception("no JIT")):
        topology = Topology._from_sbml_qual(
            model, qual, inhibited_species="S1")

    drift = topology.drift_fn(
        np.array([0.5]), np.array([1.0, 0.4]), inhibition=0.5)

    np.testing.assert_allclose(drift, np.array([-0.2]), rtol=0.0, atol=1e-12)


def test_the_roadrunner_drift_applies_scalar_inhibition() -> None:
    """The RoadRunner path scales the inhibited species by ``1 - inh``.

    The fake instance records the values assigned to it, so the test
    observes the actual modulation rather than only that ``reset`` ran.
    """
    instances: list = []

    class _FakeRoadRunner:
        def __init__(self, path):
            self.values: dict[str, float] = {}
            self.model = MagicMock()
            instances.append(self)

        def reset(self):
            self.values.clear()

        def __setitem__(self, key, value):
            self.values[key] = value

        def __getitem__(self, key):
            return self.values[key]

        def getRatesOfChange(self):
            return np.array([0.0])

    document = MagicMock()
    document.getNumErrors.return_value = 0
    model = MagicMock()
    model.getPlugin.return_value = None
    model.getListOfSpecies.return_value = [
        MagicMock(getId=lambda: "S1")]
    model.getListOfParameters.return_value = []
    document.getModel.return_value = model

    with patch.dict(
            sys.modules,
            {"roadrunner": types.SimpleNamespace(
                RoadRunner=_FakeRoadRunner)}), \
            patch("libsbml.SBMLReader.readSBML", return_value=document):
        topology = Topology.from_sbml("fake.xml", inhibited_species="S1")
        topology.drift_fn(np.array([0.5]), None, 0.5)

    assert instances[-1].values["S1"] == pytest.approx(0.5 * (1.0 - 0.5))

