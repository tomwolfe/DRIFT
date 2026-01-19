import numpy as np
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver

from concurrent.futures import ProcessPoolExecutor
import os
import cobra
import logging
from tqdm import tqdm

# Configure logging
logger = logging.getLogger(__name__)

# Global cache for workers to avoid reloading model
_worker_cache = {}


def _init_worker(model_name):
    """Initializes a worker process by loading the model once."""
    try:
        _worker_cache["solver"] = DFBASolver(model_name=model_name)
        _worker_cache["integrator"] = StochasticIntegrator(dt=0.1, noise_scale=0.03)
        _worker_cache["bridge"] = MetabolicBridge()
        logger.info(f"Worker initialized with model: {model_name}")
    except Exception as e:
        logger.error(f"Failed to initialize worker with model {model_name}: {str(e)}")
        raise


def _single_sim_wrapper(args):
    """Helper to run a single simulation in a separate process."""
    drug_kd, drug_concentration, steps, model_name = args

    try:
        # Reuse cached components if they exist
        if "solver" in _worker_cache:
            solver = _worker_cache["solver"]
            integrator = _worker_cache["integrator"]
            bridge = _worker_cache["bridge"]
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
            "time": np.arange(steps) * integrator.dt,
            "signaling": [],
            "growth": [],
            "inhibition": inhibition,
            "drug_kd": drug_kd,
        }

        for step in range(steps):
            state = integrator.step(state, inhibition)
            constraints = bridge.get_constraints(state)
            growth, _ = solver.solve_step(constraints)

            history["signaling"].append(state.copy())
            history["growth"].append(growth)

        history["signaling"] = np.array(history["signaling"])
        history["growth"] = np.array(history["growth"])
        return history
    except Exception as e:
        logger.error(f"Error in single simulation: {str(e)}")
        raise


class Workbench:
    """Multi-Scale Stochastic Research Workbench."""

    def __init__(self, drug_kd=1.0, drug_concentration=2.0, model_name="textbook"):
        """
        Initialize the Workbench with specified parameters.

        Args:
            drug_kd (float): Dissociation constant of the drug
            drug_concentration (float): Concentration of the drug in the system
            model_name (str): Name of the metabolic model to use
        """
        if drug_kd <= 0:
            raise ValueError(f"drug_kd must be positive, got {drug_kd}")
        if drug_concentration < 0:
            raise ValueError(
                f"drug_concentration must be non-negative, got {drug_concentration}"
            )

        self.binding = BindingEngine(kd=drug_kd)
        self.signaling = StochasticIntegrator(dt=0.1, noise_scale=0.03)
        self.metabolic_bridge = MetabolicBridge()
        self.solver = DFBASolver(model_name=model_name)
        self.drug_concentration = drug_concentration
        self.model_name = model_name

    def get_basal_growth(self, steps=100):
        """Calculates growth rate without any drug inhibition."""
        args = (self.binding.kd, 0.0, steps, self.model_name)
        history = _single_sim_wrapper(args)
        return np.mean(history["growth"])

    def run_simulation(self, steps=100):
        """
        Runs a single temporal simulation.

        Args:
            steps (int): Number of simulation steps to run

        Returns:
            dict: History of the simulation with time, signaling, growth, and inhibition data
        """
        if steps <= 0:
            raise ValueError(f"steps must be positive, got {steps}")

        args = (self.binding.kd, self.drug_concentration, steps, self.model_name)
        return _single_sim_wrapper(args)

    def run_monte_carlo(self, n_sims=30, steps=100, n_jobs=-1):
        """
        Runs multiple simulations with perturbed parameters in parallel.

        Args:
            n_sims (int): Number of simulations to run
            steps (int): Number of steps per simulation
            n_jobs (int): Number of parallel jobs (-1 for CPU count)

        Returns:
            dict: Results containing 'histories' and 'basal_growth'
        """
        if n_sims <= 0:
            raise ValueError(f"n_sims must be positive, got {n_sims}")
        if steps <= 0:
            raise ValueError(f"steps must be positive, got {steps}")

        # Get basal growth first for normalization
        basal_growth = self.get_basal_growth(steps=steps)

        if n_sims > 100:
            print(f"[*] Running {n_sims} Monte Carlo simulations...")

        if n_jobs == -1:
            n_jobs = os.cpu_count() or 1

        base_kd = self.binding.kd
        sim_args = []
        for _ in range(n_sims):
            perturbed_kd = base_kd * np.random.uniform(0.8, 1.2)
            sim_args.append(
                (perturbed_kd, self.drug_concentration, steps, self.model_name)
            )

        all_histories = []
        if n_jobs > 1:
            try:
                with ProcessPoolExecutor(
                    max_workers=n_jobs,
                    initializer=_init_worker,
                    initargs=(self.model_name,),
                ) as executor:
                    if n_sims > 100:
                        all_histories = list(
                            tqdm(
                                executor.map(_single_sim_wrapper, sim_args),
                                total=len(sim_args),
                                desc="Monte Carlo Simulations",
                            )
                        )
                    else:
                        all_histories = list(
                            executor.map(_single_sim_wrapper, sim_args)
                        )
            except Exception as e:
                logger.error(f"Error in multiprocessing: {str(e)}")
                n_jobs = 1  # Fallback to sequential

        if n_jobs <= 1:
            if n_sims > 100:
                for args in tqdm(sim_args, desc="Sequential Simulations"):
                    all_histories.append(_single_sim_wrapper(args))
            else:
                for args in sim_args:
                    all_histories.append(_single_sim_wrapper(args))

        return {"histories": all_histories, "basal_growth": basal_growth}
