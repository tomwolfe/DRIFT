import numpy as np
from numba import njit
from drift.topology import Topology
from drift.signaling import StochasticIntegrator
import pytest

@njit
def custom_drift_fn(state, params, feedback=None):
    """Simple drift: decay to 0.5"""
    res = 0.1 * (0.5 - state)
    return res

def test_auto_jit_custom_topology():
    # Create a topology with a jitted drift function
    topo = Topology(
        species=["S1"],
        parameters={"k": 0.1},
        drift_fn=custom_drift_fn,
        name="JitTest"
    )
    
    assert topo.jitted_step_fn is None
    
    # Initializing StochasticIntegrator should trigger auto-jit
    integrator = StochasticIntegrator(topology=topo)
    
    assert topo.jitted_step_fn is not None
    assert type(topo.jitted_step_fn).__name__ == "CPUDispatcher"
    
    # Test a step
    state = np.array([0.8])
    new_state = integrator.step(state, inhibition=0.0)
    
    assert new_state.shape == (1,)
    assert new_state[0] < 0.8 # Should decay towards 0.5

def test_no_jit_for_normal_fn():
    def normal_drift(state, params, inhibition=0.0, feedback=None):
        return 0.1 * (0.5 - state)
        
    topo = Topology(
        species=["S1"],
        parameters={"k": 0.1},
        drift_fn=normal_drift,
        name="NoJitTest"
    )
    
    integrator = StochasticIntegrator(topology=topo)
    
    # Should NOT be jitted
    assert topo.jitted_step_fn is None

if __name__ == "__main__":
    pytest.main([__file__])