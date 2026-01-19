import numpy as np
import pandas as pd
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver

from concurrent.futures import ProcessPoolExecutor
import os
import logging
from tqdm import tqdm
from typing import Dict, Any, List, Optional

# Configure logging
logger = logging.getLogger(__name__)

# Global cache for workers to avoid reloading model
_worker_cache: Dict[str, Any] = {}


def _init_worker(model_name, topology=None):
    """Initializes a worker process by loading the model once."""
    try:
        from .topology import get_default_topology
        eff_topology = topology or get_default_topology()
        _worker_cache["solver"] = DFBASolver(model_name=model_name)
        _worker_cache["integrator"] = StochasticIntegrator(dt=0.1, noise_scale=0.03, topology=eff_topology)
        _worker_cache["bridge"] = MetabolicBridge(species_names=eff_topology.species)
        logger.info(f"Worker initialized with model: {model_name}")
    except Exception as e:
        logger.error(f"Failed to initialize worker with model {model_name}: {str(e)}")
        raise


def _single_sim_wrapper(args):
    """Helper to run a single simulation in a separate process."""
    drug_kd, drug_concentration, steps, model_name, topology = args

    try:
        from .topology import get_default_topology
        eff_topology = topology or get_default_topology()

        # Reuse cached components if they exist
        if "solver" in _worker_cache:
            solver = _worker_cache["solver"]
            integrator = _worker_cache["integrator"]
            bridge = _worker_cache["bridge"]
            binding = BindingEngine(kd=drug_kd)
        else:
            # Fallback for non-pool execution
            binding = BindingEngine(kd=drug_kd)
            integrator = StochasticIntegrator(dt=0.1, noise_scale=0.03, topology=eff_topology)
            bridge = MetabolicBridge(species_names=eff_topology.species)
            solver = DFBASolver(model_name=model_name)

        inhibition = binding.calculate_inhibition(drug_concentration)
        
        # Initial state from topology
        state = integrator.topology.get_initial_state()
        feedback = 1.0  # Initial feedback (healthy state)

        history = {
            "time": np.arange(steps) * integrator.dt,
            "signaling": [],
            "growth": [],
            "inhibition": inhibition,
            "drug_kd": drug_kd,
        }

        for step in range(steps):
            # 1. Signaling step (influenced by previous metabolic feedback)
            state = integrator.step(state, inhibition, feedback=feedback)
            
            # 2. Metabolic mapping (Signaling -> Metabolism)
            constraints = bridge.get_constraints(state)
            
            # 3. FBA solver
            fba_result = solver.solve_step(constraints)
            growth = fba_result["objective_value"]
            fluxes = fba_result["fluxes"]

            # 4. Feedback mapping (Metabolism -> Signaling)
            feedback = bridge.get_feedback(fluxes)

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

    def __init__(
        self, 
        drug_kd=1.0, 
        drug_concentration=2.0, 
        model_name="textbook", 
        topology=None,
        bridge=None
    ):
        """
        Initialize the Workbench with specified parameters.

        Args:
            drug_kd (float): Dissociation constant of the drug
            drug_concentration (float): Concentration of the drug in the system
            model_name (str): Name of the metabolic model to use
            topology (Topology, optional): Signaling network topology
            bridge (MetabolicBridge, optional): Custom metabolic bridge
        """
        if drug_kd <= 0:
            raise ValueError(f"drug_kd must be positive, got {drug_kd}")
        if drug_concentration < 0:
            raise ValueError(
                f"drug_concentration must be non-negative, got {drug_concentration}"
            )

        self.binding = BindingEngine(kd=drug_kd)
        self.topology = topology
        self.signaling = StochasticIntegrator(dt=0.1, noise_scale=0.03, topology=topology)
        
        from .topology import get_default_topology
        eff_topology = topology or get_default_topology()
        
        self.metabolic_bridge = bridge or MetabolicBridge(species_names=eff_topology.species)
        # Ensure bridge has species names if it was provided without them
        if self.metabolic_bridge.species_names is None:
            self.metabolic_bridge.species_names = eff_topology.species
            
        self.solver = DFBASolver(model_name=model_name)
        self.drug_concentration = drug_concentration
        self.model_name = model_name

    def get_basal_growth(self, steps=100):
        """Calculates growth rate without any drug inhibition."""
        args = (self.binding.kd, 0.0, steps, self.model_name, self.topology)
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

        args = (self.binding.kd, self.drug_concentration, steps, self.model_name, self.topology)
        return _single_sim_wrapper(args)

    def run_monte_carlo(self, n_sims=30, steps=100, n_jobs=-1):  # noqa: C901
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

        if n_jobs == -1:
            n_jobs = os.cpu_count() or 1

        base_kd = self.binding.kd
        sim_args = []
        for _ in range(n_sims):
            perturbed_kd = base_kd * np.random.uniform(0.8, 1.2)
            sim_args.append(
                (perturbed_kd, self.drug_concentration, steps, self.model_name, self.topology)
            )

        all_histories = []
        # Progress bar logic: Always show if n_sims > 1 to provide immediate feedback
        use_tqdm = n_sims > 1

        if n_jobs > 1:
            try:
                with ProcessPoolExecutor(
                    max_workers=n_jobs,
                    initializer=_init_worker,
                    initargs=(self.model_name, self.topology),
                ) as executor:
                    if use_tqdm:
                        all_histories = list(
                            tqdm(
                                executor.map(_single_sim_wrapper, sim_args),
                                total=len(sim_args),
                                desc="Monte Carlo Simulations",
                                leave=False,
                            )
                        )
                    else:
                        all_histories = list(
                            executor.map(_single_sim_wrapper, sim_args)
                        )
            except Exception as e:
                logger.error(f"Error in multiprocessing: {str(e)}. Falling back to sequential.")
                n_jobs = 1  # Fallback to sequential

        if n_jobs <= 1:
            if use_tqdm:
                for args in tqdm(sim_args, desc="Sequential Simulations", leave=False):
                    all_histories.append(_single_sim_wrapper(args))
            else:
                for args in sim_args:
                    all_histories.append(_single_sim_wrapper(args))

        return {"histories": all_histories, "basal_growth": basal_growth}

    def export_results(self, results: Dict[str, Any], file_path: str):
        """
        Exports Monte Carlo results to a structured Parquet file.

        Args:
            results (dict): Output from run_monte_carlo
            file_path (str): Destination path (.parquet)
        """
        all_data = []
        histories = results.get("histories", [])
        species_names = self.topology.species if self.topology else ["PI3K", "AKT", "mTOR"]

        for sim_idx, hist in enumerate(histories):
            for step_idx in range(len(hist["time"])):
                row = {
                    "sim_id": sim_idx,
                    "step": step_idx,
                    "time": hist["time"][step_idx],
                    "growth": hist["growth"][step_idx],
                    "drug_kd": hist["drug_kd"],
                    "inhibition": hist.get("inhibition", 0.0),
                }
                # Add signaling species
                for i, species in enumerate(species_names):
                    row[species] = hist["signaling"][step_idx, i]
                
                all_data.append(row)

        df = pd.DataFrame(all_data)
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(file_path)), exist_ok=True)
        
        df.to_parquet(file_path, index=False)
        logger.info(f"Successfully exported {len(all_data)} rows to {file_path}")
