import numpy as np
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver

class Workbench:
    """Multi-Scale Stochastic Research Workbench."""
    def __init__(self, drug_kd=1.0, drug_concentration=2.0, model_name='textbook'):
        self.binding = BindingEngine(kd=drug_kd)
        self.signaling = StochasticIntegrator(dt=0.1, noise_scale=0.03)
        self.metabolic_bridge = MetabolicBridge()
        self.solver = DFBASolver(model_name=model_name)
        self.drug_concentration = drug_concentration

    def run_simulation(self, steps=100):
        """Runs a single temporal simulation."""
        inhibition = self.binding.calculate_inhibition(self.drug_concentration)
        # Initial state: [PI3K, AKT, mTOR]
        state = np.array([0.8, 0.8, 0.8]) 
        
        history = {
            'time': np.arange(steps) * self.signaling.dt,
            'signaling': [], # List of [PI3K, AKT, mTOR]
            'growth': [],
            'inhibition': inhibition
        }
        
        for _ in range(steps):
            state = self.signaling.step(state, inhibition)
            constraints = self.metabolic_bridge.get_constraints(state)
            growth, _ = self.solver.solve_step(constraints)
            
            history['signaling'].append(state.copy())
            history['growth'].append(growth)
            
        history['signaling'] = np.array(history['signaling'])
        history['growth'] = np.array(history['growth'])
        return history

    def run_monte_carlo(self, n_sims=30, steps=100):
        """Runs multiple simulations with perturbed parameters."""
        all_histories = []
        base_kd = self.binding.kd
        
        for _ in range(n_sims):
            # Perturb Kd by +/- 20%
            self.binding.kd = base_kd * np.random.uniform(0.8, 1.2)
            all_histories.append(self.run_simulation(steps))
            
        # Reset Kd
        self.binding.kd = base_kd
        return all_histories
