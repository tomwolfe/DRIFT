import numpy as np
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver

from concurrent.futures import ProcessPoolExecutor
import os
import cobra

# Global cache for workers
_worker_solver = None

def _init_worker(model_name):
    """Initializes a worker process by loading the model once."""
    global _worker_solver
    _worker_solver = DFBASolver(model_name=model_name)

def _single_sim_wrapper(args):
    """Helper to run a single simulation in a separate process."""
    drug_kd, drug_concentration, steps, model_name = args
    # Use the cached solver if available, otherwise create one (fallback)
    global _worker_solver
    if _worker_solver is None:
        wb = Workbench(drug_kd=drug_kd, drug_concentration=drug_concentration, model_name=model_name)
    else:
        # Create a light workbench that reuses the process's solver
        wb = Workbench(drug_kd=drug_kd, drug_concentration=drug_concentration, model_name=model_name)
        wb.solver = _worker_solver
        
    return wb.run_simulation(steps)

class Workbench:
    """Multi-Scale Stochastic Research Workbench."""
    def __init__(self, drug_kd=1.0, drug_concentration=2.0, model_name='textbook'):
        self.binding = BindingEngine(kd=drug_kd)
        self.signaling = StochasticIntegrator(dt=0.1, noise_scale=0.03)
        self.metabolic_bridge = MetabolicBridge()
        self.solver = DFBASolver(model_name=model_name)
        self.drug_concentration = drug_concentration
        self.model_name = model_name

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

    def run_monte_carlo(self, n_sims=30, steps=100, n_jobs=-1):
        """Runs multiple simulations with perturbed parameters in parallel."""
        if n_jobs == -1:
            n_jobs = os.cpu_count() or 1
            
        base_kd = self.binding.kd
        sim_args = []
        
        for _ in range(n_sims):
            perturbed_kd = base_kd * np.random.uniform(0.8, 1.2)
            sim_args.append((perturbed_kd, self.drug_concentration, steps, self.model_name))
            
        if n_jobs > 1:
            with ProcessPoolExecutor(max_workers=n_jobs, initializer=_init_worker, initargs=(self.model_name,)) as executor:
                all_histories = list(executor.map(_single_sim_wrapper, sim_args))
        else:
            all_histories = []
            for args in sim_args:
                all_histories.append(_single_sim_wrapper(args))
                
        return all_histories
