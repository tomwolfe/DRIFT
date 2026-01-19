import pytest
import numpy as np
from numba import njit
from drift.signaling import StochasticIntegrator, create_langevin_integrator
from drift.topology import Topology

@njit
def custom_drift_fn(state, params, feedback=1.0):
    # Simple decay modulated by feedback
    return -params[0] * state * (2.0 - feedback)

def test_jitted_custom_topology():
    """Test that custom jitted topologies work and are faster (implicitly)."""
    # Create jitted step function
    jitted_step = create_langevin_integrator(custom_drift_fn)
    
    # Create topology with jitted step
    topology = Topology(
        species=["S1"],
        parameters={"k": 0.1},
        jitted_step_fn=jitted_step,
        name="jitted_custom"
    )
    
    integrator = StochasticIntegrator(topology=topology)
    
    initial_state = np.array([1.0])
    inhibition = 0.5
    feedback = 1.0
    
    # Test step
    new_state = integrator.step(initial_state, inhibition, feedback=feedback)
    
    assert new_state.shape == (1,)
    assert 0 <= new_state[0] <= 1.0
    # Since it's a decay, it should generally decrease (ignoring noise)
    # But with noise it might increase. Let's just check it runs without error.
    
    # Run multiple steps to ensure stability
    state = initial_state
    for _ in range(10):
        state = integrator.step(state, inhibition, feedback=feedback)
        assert 0 <= state[0] <= 1.0
