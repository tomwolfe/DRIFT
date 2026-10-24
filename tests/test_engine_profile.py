"""The PK seam: DRIFT's engine driven by a concentration-TIME PROFILE.

The two engines are ONE trial only if the plasma concentration profile
is what the PD layer sees.  These tests pin that seam:

  * a CONSTANT profile is the same run as the scalar concentration --
    the profile path is not a different model, it is the same model with
    the exposure read at the step's own time;
  * the engine RECORDS the exposure it was given: the concentration
    series is the profile sampled at the step times, and the inhibition
    series is the Hill occupancy of that profile;
  * a profile that falls and rises produces an inhibition series that
    falls and rises -- the exposure is read per step, not summarised
    once;
  * a malformed profile is REFUSED rather than silently reinterpreted.

Every property is MEASURED from the engine's own history.
"""

from __future__ import annotations

import numpy as np
import pytest

from drift.binding import BindingEngine
from drift.engine import SimulationEngine
from drift.metabolic import DFBASolver, MetabolicBridge
from drift.signaling import StochasticIntegrator
from drift.topology import get_default_topology

DT = 0.05


def _make_engine(kd: float = 0.5, hill: float = 1.5) -> SimulationEngine:
    integrator = StochasticIntegrator(
        dt=DT, noise_scale=0.01, topology=get_default_topology())
    solver = DFBASolver(model_name="textbook")
    binding = BindingEngine(targets=kd, hill_coefficient=hill)
    return SimulationEngine(integrator, solver, MetabolicBridge(), binding)


def _hill(concentration: float, kd: float, hill: float) -> float:
    if concentration <= 0.0:
        return 0.0
    numerator = concentration ** hill
    return float(numerator / (numerator + kd ** hill))


def test_a_constant_profile_is_the_scalar_run() -> None:
    """The profile path is the scalar path with a constant exposure."""
    concentration = 2.0
    profile = np.full(20, concentration, dtype=np.float64)
    times = np.arange(20, dtype=np.float64) * DT

    scalar = _make_engine().run(steps=20, drug_concentration=concentration, seed=11)
    driven = _make_engine().run(steps=20, drug_profile=profile,
                                drug_profile_times=times, seed=11)

    assert np.allclose(driven["growth"], scalar["growth"], rtol=0.0, atol=0.0)
    assert driven["cell_death"] == scalar["cell_death"]
    assert driven["inhibition"] == pytest.approx(scalar["inhibition"])


def test_the_engine_records_the_profile_it_was_given() -> None:
    """The engine's own record IS the plasma profile and its occupancy."""
    kd, hill = 0.5, 1.5
    profile = np.array([0.0, 0.5, 2.0, 6.0, 3.0, 1.0], dtype=np.float64)
    times = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], dtype=np.float64)
    steps = 10

    history = _make_engine(kd, hill).run(
        steps=steps, drug_profile=profile, drug_profile_times=times, seed=5)

    expected_concentration = np.interp(
        np.arange(steps) * DT, times, profile)
    assert np.allclose(history["concentration"], expected_concentration,
                       rtol=0.0, atol=1e-12)
    expected_inhibition = np.asarray(
        [_hill(float(c), kd, hill) for c in expected_concentration])
    assert np.allclose(history["inhibition_series"], expected_inhibition,
                       rtol=0.0, atol=1e-12)
    # The reported inhibition is the PEAK of the profile's exposure.
    assert history["inhibition"] == pytest.approx(
        float(expected_inhibition.max()))


def test_the_exposure_is_read_at_every_step() -> None:
    """A profile that falls and rises gives an inhibition that falls and
    rises -- the exposure is not summarised once."""
    kd, hill = 0.5, 1.5
    profile = np.array([4.0, 1.0, 4.0], dtype=np.float64)
    times = np.array([0.0, 0.05, 0.10], dtype=np.float64)

    history = _make_engine(kd, hill).run(
        steps=3, drug_profile=profile, drug_profile_times=times, seed=5)

    occupancy = _hill(4.0, kd, hill)
    assert history["inhibition_series"][0] == pytest.approx(occupancy)
    assert history["inhibition_series"][1] == pytest.approx(
        _hill(1.0, kd, hill))
    assert history["inhibition_series"][2] == pytest.approx(occupancy)
    assert history["inhibition_series"][1] < history["inhibition_series"][0]


def test_the_profile_may_be_given_without_times() -> None:
    """Without times the profile is read on the engine's own time grid."""
    profile = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    history = _make_engine().run(steps=3, drug_profile=profile, seed=5)
    assert np.allclose(history["concentration"], profile, rtol=0.0, atol=0e-12)


def test_a_malformed_profile_is_refused() -> None:
    """A malformed exposure is refused, not silently reinterpreted."""
    engine = _make_engine()
    with pytest.raises(ValueError):
        engine.run(steps=3, drug_profile=np.array([1.0, 2.0]),
                   drug_profile_times=np.array([0.0]))
    with pytest.raises(ValueError):
        engine.run(steps=3, drug_profile=np.array([1.0, 2.0]),
                   drug_profile_times=np.array([0.5, 0.0]))
    with pytest.raises(ValueError):
        engine.run(steps=3, drug_profile=np.array([1.0, -2.0]))
    with pytest.raises(ValueError):
        engine.run(steps=3, drug_profile=np.array([]))
    with pytest.raises(ValueError):
        engine.run(steps=3, drug_profile_times=np.array([0.0]))


def test_a_run_without_any_exposure_is_refused() -> None:
    """Neither a scalar nor a profile is not a run."""
    with pytest.raises(ValueError):
        _make_engine().run(steps=3)
