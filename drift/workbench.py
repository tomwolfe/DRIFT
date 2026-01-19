import numpy as np
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver

from concurrent.futures import ProcessPoolExecutor
import os
import cobra

# Global cache for workers to avoid reloading model
_worker_cache = {}

def _init_worker(model_name):
    """Initializes a worker process by loading the model once."""
    global _worker_cache
    _worker_cache['solver'] = DFBASolver(model_name=model_name)
    _worker_cache['integrator'] = StochasticIntegrator(dt=0.1, noise_scale=0.03)
    _worker_cache['bridge'] = MetabolicBridge()

def _single_sim_wrapper(args):
    """Helper to run a single simulation in a separate process."""
    drug_kd, drug_concentration, steps, model_name = args
    global _worker_cache
    
    # Reuse cached components if they exist
    if 'solver' in _worker_cache:
        solver = _worker_cache['solver']
        integrator = _worker_cache['integrator']
        bridge = _worker_cache['bridge']
        binding = BindingEngine(kd=drug_kd)
    else:
        # Fallback for non-pool execution
        binding = BindingEngine(kd=drug_kd)
        integrator = StochasticIntegrator(dt=0.1, noise_scale=0.03)
        bridge = MetabolicBridge()
        solver = DFBASolver(model_name=model_name)
        
    inhibition = binding.calculate_inhibition(drug_concentration)
    state = np.array([0.8, 0.8, 0.8]) 
    
    history = {
        'time': np.arange(steps) * integrator.dt,
        'signaling': [],
        'growth': [],
        'inhibition': inhibition
    }
    
    for _ in range(steps):
        state = integrator.step(state, inhibition)
        constraints = bridge.get_constraints(state)
        growth, _ = solver.solve_step(constraints)
        
        history['signaling'].append(state.copy())
        history['growth'].append(growth)
        
    history['signaling'] = np.array(history['signaling'])
    history['growth'] = np.array(history['growth'])
    return history

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
        # This now just calls the same logic as the wrapper for consistency
        args = (self.binding.kd, self.drug_concentration, steps, self.model_name)
        return _single_sim_wrapper(args)

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
