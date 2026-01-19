import numpy as np
import pytest
from drift.signaling import StochasticIntegrator
from drift.topology import Topology, drift_model
from scipy.integrate import odeint

def test_stochastic_mean_convergence():
    """
    Scientific Validation: Verify that the stochastic mean of the Milstein scheme
    converges to the deterministic ODE solution as noise scale decreases or N increases.
    """
    # 1. Define a simple linear decay model
    # dx = (k_base - k_deg * x) * dt + noise * sqrt(x * (1-x)) * dW
    
    @drift_model("linear_decay")
    def linear_drift(state, params, feedback=None):
        k_base = params[0]
        k_deg = params[1]
        inhibition = params[2]
        
        # Simple linear dynamics: dx/dt = k_base - k_deg * x
        # Apply inhibition to the production rate
        return np.array([k_base * (1.0 - inhibition) - k_deg * state[0]])

    topo = Topology(
        species=["X"],
        parameters={"k_base": 0.2, "k_deg": 0.2},
        drift_fn=linear_drift,
        name="LinearDecay"
    )

    # 2. Deterministic Solution using scipy.integrate.odeint
    def ode_func(x, t, k_base, k_deg, inhibition):
        return k_base * (1.0 - inhibition) - k_deg * x

    t_eval = np.linspace(0, 10, 101)
    k_base, k_deg = 0.2, 0.2
    inhibition = 0.5
    x0 = [0.5]
    
    sol_ode = odeint(ode_func, x0[0], t_eval, args=(k_base, k_deg, inhibition))
    
    # 3. Stochastic Ensemble
    # We use a small noise scale and large N to see convergence to mean
    dt = 0.1
    noise_scale = 0.01 
    integrator = StochasticIntegrator(dt=dt, noise_scale=noise_scale, topology=topo)
    
    n_sims = 100
    steps = 100 # 10 seconds total (100 * 0.1)
    
    ensemble_results = []
    for _ in range(n_sims):
        state = np.array([0.5])
        history = [state[0]]
        for _ in range(steps):
            state = integrator.step(state, inhibition=inhibition)
            history.append(state[0])
        ensemble_results.append(history)
    
    stochastic_mean = np.mean(ensemble_results, axis=0)
    
    # 4. Compare
    # The mean of the stochastic process should be very close to the ODE solution
    # Especially for linear systems where E[f(x)] = f(E[x])
    mse = np.mean((sol_ode.flatten() - stochastic_mean)**2)
    
    print(f"Mean Squared Error between ODE and Stochastic Mean: {mse:.6f}")
    
    # Allow for some statistical fluctuation, but it should be small
    assert mse < 0.005, f"Stochastic mean deviated too much from ODE (MSE={mse})"

if __name__ == "__main__":
    test_stochastic_mean_convergence()
