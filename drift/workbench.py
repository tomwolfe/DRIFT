import numpy as np
import pandas as pd
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver
from .engine import SimulationEngine, SimulationResult

from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import os
import logging
import json
import random
from tqdm import tqdm
from typing import Dict, Any, List, Optional, Union, Generator

# Configure logging
logger = logging.getLogger(__name__)

# Global cache for workers to avoid reloading model
# Key format: (model_fingerprint, topology_name, bridge_hash)
_worker_cache: Dict[tuple[str, str, int], Any] = {}


def clear_worker_cache():
    """Clears the global worker cache. Useful for testing or when models change."""
    _worker_cache.clear()
    logger.info("Worker cache cleared.")


def _get_bridge_hash(bridge: Optional[MetabolicBridge]) -> int:
    """Generates a stable hash for a MetabolicBridge configuration."""
    if bridge is None:
        return 0
    
    # Extract structural data for hashing
    try:
        mapping_data = []
        if hasattr(bridge, "mappings"):
            for m in bridge.mappings:
                # Handle callable mapping_fn by using its name or code hash
                m_copy = {k: (str(v) if callable(v) else v) for k, v in m.items()}
                mapping_data.append(m_copy)
        
        rev_mapping_data = []
        if hasattr(bridge, "reverse_mappings"):
            for rm in bridge.reverse_mappings:
                rm_copy = {k: (str(v) if callable(v) else v) for k, v in rm.items()}
                rev_mapping_data.append(rm_copy)
        
        # Include other bridge state that affects results
        state_data = {
            "basal_growth_rate": getattr(bridge, "basal_growth_rate", 0.2),
            "strict_mapping": getattr(bridge, "strict_mapping", False),
            "species_names": getattr(bridge, "species_names", [])
        }
        
        # Use a deterministic string representation for hashing
        hash_str = json.dumps([mapping_data, rev_mapping_data, state_data], sort_keys=True)
        import hashlib
        return int(hashlib.sha256(hash_str.encode()).hexdigest(), 16) % (2**32)  # Limit to 32-bit integer
    except Exception as e:
        logger.warning(f"Failed to generate deterministic hash for bridge: {e}. Falling back to default hash.")
        return hash(str(getattr(bridge, "__dict__", "")))


def _init_worker(model_or_name, topology=None, bridge=None, model_fingerprint=None):
    """Initializes a worker process by loading the model once."""
    try:
        from .topology import get_default_topology
        eff_topology = topology or get_default_topology()
        
        # If fingerprint not provided, try to get it
        if model_fingerprint is None:
            if isinstance(model_or_name, str):
                # Temporary solver to get fingerprint
                temp_solver = DFBASolver(model_name=model_or_name)
                model_fingerprint = temp_solver.get_fingerprint()
            else:
                # It's already a model object
                temp_solver = DFBASolver(model=model_or_name)
                model_fingerprint = temp_solver.get_fingerprint()

        # Create a unique key for this configuration
        bridge_hash = _get_bridge_hash(bridge)
        cache_key = (model_fingerprint, eff_topology.name, bridge_hash)
        
        if cache_key not in _worker_cache:
            if isinstance(model_or_name, str):
                solver = DFBASolver(model_name=model_or_name)
            else:
                solver = DFBASolver(model=model_or_name)

            _worker_cache[cache_key] = {
                "solver": solver,
                "integrator": StochasticIntegrator(dt=0.1, noise_scale=0.03, topology=eff_topology),
                "bridge": bridge or MetabolicBridge(species_names=eff_topology.species)
            }
            logger.info(f"Worker {os.getpid()} initialized with model fingerprint: {model_fingerprint[:8]}..., topology: {eff_topology.name}")
    except Exception as e:
        logger.error(f"Failed to initialize worker: {str(e)}")
        raise


def _single_sim_wrapper(args) -> SimulationResult:
    """Helper to run a single simulation in a separate process."""
    drug_kd, drug_concentration, steps, model_or_name, topology, bridge, custom_params, sub_steps, refine_feedback, model_fingerprint, seed = args

    try:
        from .topology import get_default_topology
        eff_topology = topology or get_default_topology()
        bridge_hash = _get_bridge_hash(bridge)
        
        # If fingerprint not provided, we have to calculate it (expensive in worker)
        if model_fingerprint is None:
             if isinstance(model_or_name, str):
                temp_solver = DFBASolver(model_name=model_or_name)
                model_fingerprint = temp_solver.get_fingerprint()
             else:
                temp_solver = DFBASolver(model=model_or_name)
                model_fingerprint = temp_solver.get_fingerprint()

        cache_key = (model_fingerprint, eff_topology.name, bridge_hash)

        # Reuse cached components if they exist
        if cache_key in _worker_cache:
            cached = _worker_cache[cache_key]
            solver = cached["solver"]
            integrator = cached["integrator"]
            sim_bridge = bridge or cached["bridge"]
            binding = BindingEngine(targets=drug_kd)
        else:
            # Fallback for non-pool execution or if init_worker missed this key
            binding = BindingEngine(targets=drug_kd)
            integrator = StochasticIntegrator(dt=0.1, noise_scale=0.03, topology=eff_topology)
            sim_bridge = bridge or MetabolicBridge(species_names=eff_topology.species)
            
            if isinstance(model_or_name, str):
                solver = DFBASolver(model_name=model_or_name)
            else:
                solver = DFBASolver(model=model_or_name)

        engine = SimulationEngine(
            integrator=integrator,
            solver=solver,
            bridge=sim_bridge,
            binding=binding
        )

        return engine.run(
            steps=steps,
            drug_concentration=drug_concentration,
            sub_steps=sub_steps,
            refine_feedback=refine_feedback,
            custom_params=custom_params,
            seed=seed
        )
    except Exception as e:
        logger.error(f"Error in single simulation: {str(e)}")
        raise


class Workbench:
    """Multi-Scale Stochastic Research Workbench."""

    def __init__(
        self, 
        drug_kd: Union[float, Dict[str, float]] = 1.0, 
        drug_concentration: float = 2.0, 
        model_name="textbook", 
        topology=None,
        bridge=None,
        time_unit="hours",
        concentration_unit="uM",
        flux_unit="mmol/gDW/h",
        sub_steps=5,
        refine_feedback=False,
        strict_mapping=False,
        random_state: Optional[int] = None
    ):
        """
        Initialize the Workbench with specified parameters.

        Args:
            drug_kd (float or dict): Dissociation constant(s) of the drug (in concentration_units)
            drug_concentration (float): Concentration of the drug in the system (in concentration_units)
            model_name (str): Name of the metabolic model to use
            topology (Topology, optional): Signaling network topology
            bridge (MetabolicBridge, optional): Custom metabolic bridge
            time_unit (str): Unit for time (default: "hours")
            concentration_unit (str): Unit for concentration (default: "uM")
            flux_unit (str): Unit for metabolic flux (default: "mmol/gDW/h")
            sub_steps (int): Number of signaling sub-steps per metabolic step (default: 5)
            refine_feedback (bool): If True, uses a predictor-corrector step for metabolic feedback.
            strict_mapping (bool): If True, enforces strict reaction ID mapping in the bridge.
            random_state (int, optional): Seed for reproducibility.
        """
        # Validate drug_kd
        if isinstance(drug_kd, (int, float)):
            if drug_kd <= 0:
                raise ValueError(f"drug_kd must be positive, got {drug_kd}")
        elif isinstance(drug_kd, dict):
            for name, val in drug_kd.items():
                if val <= 0:
                    raise ValueError(f"drug_kd for {name} must be positive, got {val}")
        else:
            raise TypeError("drug_kd must be float or dict")

        if drug_concentration < 0:
            raise ValueError(
                f"drug_concentration must be non-negative, got {drug_concentration}"
            )

        self.random_state = random_state
        if random_state is not None:
            np.random.seed(random_state)
            random.seed(random_state)

        self.binding = BindingEngine(targets=drug_kd)
        self.topology = topology
        self.signaling = StochasticIntegrator(dt=0.1, noise_scale=0.03, topology=topology)
        
        from .topology import get_default_topology
        eff_topology = topology or get_default_topology()
        
        self.metabolic_bridge = bridge or MetabolicBridge(
            species_names=eff_topology.species,
            strict_mapping=strict_mapping,
            flux_unit=flux_unit
        )
        if self.metabolic_bridge.species_names is None:
            self.metabolic_bridge.species_names = eff_topology.species
            
        self.solver = DFBASolver(model_name=model_name)
        if self.solver.headless:
            print("\n" + "!"*80)
            print("!!! WARNING: DFBASolver is in HEADLESS mode. No COBRA-compatible LP solver found. !!!")
            print("!!! Metabolism is PROXIED (QUALITATIVE ONLY). Signaling will run with PROXIED feedback. !!!")
            print("!!! To fix this, install a solver: pip install swiglpk                          !!!")
            print("!"*80 + "\n")
            logger.warning("Workbench initialized with headless DFBASolver (QUALITATIVE ONLY).")
        else:
            # Validate bridge against model IDs
            self.metabolic_bridge.validate_with_model(self.solver.model)
            
            # Pareto: Auto-calibrate normalization if we have a model
            try:
                logger.info("Auto-calibrating MetabolicBridge basal growth rate...")
                if self.solver.model is not None:
                    basal_sol = self.solver.model.optimize()
                    if basal_sol is not None and basal_sol.status == "optimal":
                        self.metabolic_bridge.calibrate(basal_sol.fluxes.to_dict())
            except Exception as e:
                logger.warning(f"Auto-calibration failed: {e}")

        self.drug_concentration = drug_concentration
        self.model_name = model_name
        self.time_unit = time_unit
        self.concentration_unit = concentration_unit
        self.flux_unit = flux_unit
        self.sub_steps = sub_steps
        self.refine_feedback = refine_feedback

    def get_basal_growth(self, steps=100) -> float:
        """Calculates growth rate without any drug inhibition."""
        fingerprint = self.solver.get_fingerprint()
        seed = self.random_state if self.random_state is not None else None
        args = (self.binding.kd, 0.0, steps, self.solver.model or self.model_name, self.topology, self.metabolic_bridge, None, self.sub_steps, self.refine_feedback, fingerprint, seed)
        history = _single_sim_wrapper(args)
        if isinstance(history, dict) and "growth" in history:
            return float(np.mean(history["growth"]))
        else:
            logger.warning("History is not a dictionary with 'growth' key, returning 0.0")
            return 0.0

    def run_simulation(self, steps=100, seed: Optional[int] = None) -> SimulationResult:
        """
        Runs a single temporal simulation.

        Args:
            steps (int): Number of simulation steps to run
            seed (int, optional): Seed for this specific simulation.

        Returns:
            SimulationResult: History of the simulation
        """
        if steps <= 0:
            raise ValueError(f"steps must be positive, got {steps}")

        eff_seed = seed if seed is not None else self.random_state
        fingerprint = self.solver.get_fingerprint()
        args = (self.binding.kd, self.drug_concentration, steps, self.solver.model or self.model_name, self.topology, self.metabolic_bridge, None, self.sub_steps, self.refine_feedback, fingerprint, eff_seed)
        return _single_sim_wrapper(args)

    def _history_to_rows(self, hist, sim_idx):
        """Internal helper to convert a history dict to a list of flat rows."""
        rows = []
        from .topology import get_default_topology
        eff_topology = self.topology or get_default_topology()
        species_names = eff_topology.species

        for step_idx in range(len(hist["time"])):
            row = {
                "sim_id": sim_idx,
                "step": step_idx,
                "time": hist["time"][step_idx],
                "time_unit": self.time_unit,
                "growth": hist["growth"][step_idx],
                "drug_kd": hist["drug_kd"],
                "concentration_unit": self.concentration_unit,
                "inhibition": hist.get("inhibition", 0.0),
            }
            # Add signaling species
            for i, species in enumerate(species_names):
                row[species] = hist["signaling"][step_idx, i]
            
            # Add custom parameters
            for p_name, p_val in hist.get("params", {}).items():
                row[f"param_{p_name}"] = p_val
            
            rows.append(row)
        return rows

    def run_monte_carlo(self, n_sims=30, steps=100, n_jobs=-1, perturb_params: Optional[List[str]] = None, export_to: Optional[str] = None, return_generator: bool = False) -> Union[Dict[str, Any], Generator[SimulationResult, None, None]]:  # noqa: C901
        """
        Runs multiple simulations with perturbed parameters in parallel.

        Args:
            n_sims (int): Number of simulations to run
            steps (int): Number of steps per simulation
            n_jobs (int): Number of parallel jobs (-1 for CPU count)
            perturb_params (Optional[List[str]]): List of signaling parameter names to perturb.
            export_to (Optional[str], optional): Path to a Parquet file for incremental writing.
            return_generator (bool): If True, returns a generator yielding histories one by one.

        Returns:
            dict or generator: Results containing 'histories' or a generator of histories.
        """
        if n_sims <= 0:
            raise ValueError(f"n_sims must be positive, got {n_sims}")
        if steps <= 0:
            raise ValueError(f"steps must be positive, got {steps}")

        # Memory Safety Enforcement (Pareto: Prevent OOM by default for huge ensembles)
        total_data_points = n_sims * steps
        if total_data_points > 2_000_000 and not (export_to or return_generator):
            logger.warning(
                f"Large ensemble ({total_data_points} pts) detected. "
                "Automatically switching to generator mode to prevent OOM."
            )
            return_generator = True
        elif total_data_points > 100_000 and not (export_to or return_generator):
            logger.info(
                f"Moderate simulation ensemble ({total_data_points} data points) "
                "monitored for memory pressure."
            )

        print(f"[*] Establishing basal phenotypic baseline...")
        basal_growth = self.get_basal_growth(steps=steps)
        print(f"[+] Basal Growth: {basal_growth:.4f} {self.time_unit}⁻¹")

        if n_jobs == -1:
            n_jobs = os.cpu_count() or 1

        base_kd = self.binding.targets
        sim_args = []
        from .topology import get_default_topology
        eff_topology = self.topology or get_default_topology()
        
        fingerprint = self.solver.get_fingerprint()
        model_or_name = self.solver.model or self.model_name

        # Master seeding logic
        master_seed = self.random_state if self.random_state is not None else np.random.randint(0, 2**31)
        sim_seeds = [master_seed + i for i in range(n_sims)]

        for i in range(n_sims):
            # Ensure parameter perturbation is also deterministic if random_state is set
            if self.random_state is not None:
                np.random.seed(sim_seeds[i])

            if isinstance(base_kd, dict):
                perturbed_kd = {name: kd * np.random.uniform(0.8, 1.2) for name, kd in base_kd.items()}
            else:
                # Should not happen with new BindingEngine but for safety
                perturbed_kd = float(base_kd) * np.random.uniform(0.8, 1.2)

            custom_params = {}
            if perturb_params:
                for p_name in perturb_params:
                    if p_name in eff_topology.parameters:
                        custom_params[p_name] = eff_topology.parameters[p_name] * np.random.uniform(0.8, 1.2)

            sim_args.append(
                (perturbed_kd, self.drug_concentration, steps, model_or_name, self.topology, self.metabolic_bridge, custom_params, self.sub_steps, self.refine_feedback, fingerprint, sim_seeds[i])
            )

        print(f"[*] Starting Monte Carlo Ensemble (N={n_sims}, Workers={n_jobs})...")

        def _gen_results():
            writer = None
            try:
                if n_jobs > 1:
                    # 'spawn' is safer than 'fork' for processes with complex C-extensions (like COBRA solvers)
                    ctx = mp.get_context('spawn')
                    with ProcessPoolExecutor(
                        max_workers=n_jobs,
                        mp_context=ctx,
                        initializer=_init_worker,
                        initargs=(model_or_name, self.topology, self.metabolic_bridge, fingerprint),
                    ) as executor:
                        iterator = executor.map(_single_sim_wrapper, sim_args)
                        for i, hist in enumerate(tqdm(iterator, total=n_sims, desc="Monte Carlo", disable=n_sims < 2)):
                            if export_to:
                                df_chunk = pd.DataFrame(self._history_to_rows(hist, i))
                                if writer is None:
                                    import pyarrow as pa
                                    import pyarrow.parquet as pq
                                    table = pa.Table.from_pandas(df_chunk)
                                    writer = pq.ParquetWriter(export_to, table.schema)
                                    writer.write_table(table)
                                else:
                                    table = pa.Table.from_pandas(df_chunk)
                                    writer.write_table(table)
                            yield hist
                else:
                    for i, args in enumerate(tqdm(sim_args, desc="Sequential MC", disable=n_sims < 2)):
                        hist = _single_sim_wrapper(args)
                        if export_to:
                            df_chunk = pd.DataFrame(self._history_to_rows(hist, i))
                            if writer is None:
                                import pyarrow as pa
                                import pyarrow.parquet as pq
                                table = pa.Table.from_pandas(df_chunk)
                                writer = pq.ParquetWriter(export_to, table.schema)
                                writer.write_table(table)
                            else:
                                table = pa.Table.from_pandas(df_chunk)
                                writer.write_table(table)
                        yield hist
            finally:
                if writer:
                    writer.close()
                    print(f"[+] Incremental export saved to {export_to}")

        if return_generator or export_to:
            return _gen_results()  # type: ignore

        # Default behavior: collect all in memory
        # SAFETY CHECK: If the user is about to crash their machine, force generator
        if n_sims > 1000:
             logger.warning(f"Extremely large N={n_sims} requested without export_to. Forcing generator mode to prevent OOM.")
             return _gen_results()  # type: ignore

        all_histories = list(_gen_results())
        
        dead_count = sum(1 for h in all_histories if h.get("cell_death", False))
        if dead_count > 0:
            print(f"[!] Warning: {dead_count}/{n_sims} simulations resulted in metabolic collapse (Cell Death).")

        return {"histories": all_histories, "basal_growth": basal_growth}

    def run_global_sensitivity_analysis(self, n_sims=100, steps=100, n_jobs=-1):
        """
        Performs Global Sensitivity Analysis (GSA) to identify drivers of metabolic drift.
        Uses both Spearman correlation and Sobol-inspired Variance Decomposition (S1).
        """
        print(f"[*] Starting Global Sensitivity Analysis (N={n_sims})...")
        
        # Identify all parameters to perturb
        from .topology import get_default_topology
        eff_topology = self.topology or get_default_topology()
        params_to_perturb = list(eff_topology.parameters.keys())
        
        # We need more simulations for Sobol than for simple correlation
        if n_sims < len(params_to_perturb) * 2:
            logger.warning(f"N={n_sims} may be too low for robust Sobol indices. Recommended: >{len(params_to_perturb)*10}")

        results = self.run_monte_carlo(n_sims=n_sims, steps=steps, n_jobs=n_jobs, perturb_params=params_to_perturb)

        # 1. Collect Data
        data = []
        if isinstance(results, dict) and "histories" in results:
            histories = results["histories"]
        else:
            histories = list(results) # Generator

        for h in histories:
            if isinstance(h, dict):
                # Calculate aggregate metrics for sensitivity
                avg_growth = np.mean(h["growth"])
                # We also look at "drift" (instability)
                growth_cv = np.std(h["growth"]) / (avg_growth + 1e-6)
                
                row = {
                    "growth": avg_growth, 
                    "instability": growth_cv,
                    "drug_kd": h["drug_kd"]
                }
                row.update(h["params"])
                data.append(row)
        
        df = pd.DataFrame(data)
        
        # 2. Linear/Monotonic Sensitivity (Spearman)
        correlation_matrix = df.corr(method="spearman")
        growth_corr = correlation_matrix["growth"].drop("growth", errors="ignore")
        correlations = growth_corr.sort_values(ascending=False)

        print("\n[+] Monotonic Sensitivity (Spearman Correlation with Growth):")
        for param, corr in correlations.items():
            if param in params_to_perturb or param == "drug_kd":
                print(f"    - {param:20}: {corr: .4f}")

        # 3. Variance-based Sensitivity (Sobol S1)
        # S1 = Var(E[Y|Xi]) / Var(Y)
        sobol_s1 = {}
        total_var = df["growth"].var()
        
        if total_var > 1e-9:
            for param in params_to_perturb + ["drug_kd"]:
                if param in df.columns:
                    # Use binning to estimate E[Y|Xi]
                    # Pareto: 10 bins or sqrt(N) for better estimation
                    n_bins = int(np.sqrt(len(df)))
                    if n_bins < 2: n_bins = 2
                    if n_bins > 10: n_bins = 10
                    
                    try:
                        df["bin"] = pd.qcut(df[param], n_bins, labels=False, duplicates='drop')
                        conditional_means = df.groupby("bin")["growth"].mean()
                        if len(conditional_means) > 1:
                            var_of_means = conditional_means.var()
                            sobol_s1[param] = var_of_means / total_var
                    except Exception as e:
                        logger.debug(f"Sobol estimation failed for {param}: {e}")

            print("\n[+] Variance-based Sensitivity (Sobol S1 Indices):")
            sorted_s1 = sorted(sobol_s1.items(), key=lambda x: x[1], reverse=True)
            for param, s1 in sorted_s1:
                print(f"    - {param:20}: {s1:.4f}")

        # 4. Interaction Detection (Total Effect vs First Order)
        # If sum(S1) < 1.0 significantly, it suggests strong interactions (S2, S3...)
        sum_s1 = sum(sobol_s1.values())
        if sum_s1 < 0.7:
             print(f"\n[!] High Interaction Detected: Sum(S1) = {sum_s1:.2f}")
             print("    Metabolic drift is driven by non-linear parameter couplings.")

        return {
            "spearman": correlations.to_dict(),
            "sobol_s1": sobol_s1,
            "sum_s1": sum_s1
        }

    def calibrate_to_data(
        self, 
        experimental_data_path: str, 
        parameter_name: str = "drug_kd", 
        param_range: tuple = (0.01, 10.0), 
        n_points: int = 10,
        mapping: Optional[Dict[str, str]] = None
    ):
        """
        Calibrates a specific parameter by minimizing MSE against experimental data.
        
        Args:
            experimental_data_path: Path to experimental CSV.
            parameter_name: Name of parameter to calibrate (e.g., 'drug_kd').
            param_range: (min, max) for the parameter.
            n_points: Number of points to sample in the range (Grid Search).
            mapping: Mapping from experimental columns to simulation keys.
        """
        from .validation import ExternalValidator
        validator = ExternalValidator(experimental_data_path)
        
        if mapping is None:
            mapping = {"growth": "growth", "pAKT": "AKT"}

        print(f"[*] Calibrating {parameter_name} against {experimental_data_path}...")
        
        search_space = np.linspace(param_range[0], param_range[1], n_points)
        best_param = None
        min_total_mse = float("inf")
        
        for val in search_space:
            # Update parameter
            if parameter_name == "drug_kd":
                self.binding.kd = val
            else:
                # Handle topology parameters if needed
                pass
                
            # Run simulation
            # Experimental data often covers long time scales, match steps to data
            if validator.exp_data is not None:
                max_time = validator.exp_data["time"].max()
                steps = int(max_time / self.signaling.dt)
            else:
                logger.error("Experimental data not loaded for calibration.")
                raise ValueError("Experimental data not loaded for calibration.")
            
            history = self.run_simulation(steps=steps)
            
            # Benchmark
            try:
                # Cast SimulationResult to dict to satisfy mypy
                history_dict = dict(history)
                metrics = validator.benchmark(history_dict, mapping)
                total_mse = sum(m["mse"] for m in metrics.values())
                
                if total_mse < min_total_mse:
                    min_total_mse = total_mse
                    best_param = val
                    
                logger.info(f"  - {parameter_name}={val:.4f} | Total MSE: {total_mse:.6f}")
            except Exception as e:
                logger.warning(f"  - {parameter_name}={val:.4f} | Calibration step failed: {e}")

        if best_param is not None:
            print(f"[+] Calibration Complete. Best {parameter_name}: {best_param:.4f} (MSE: {min_total_mse:.6f})")
            if parameter_name == "drug_kd":
                self.binding.kd = best_param
        
        return best_param, min_total_mse

    def export_results(self, results: Dict[str, Any], file_path: str):
        """
        Exports Monte Carlo results to a structured Parquet file.
        """
        all_data = []
        histories = results.get("histories", [])
        
        for sim_idx, hist in enumerate(histories):
            all_data.extend(self._history_to_rows(hist, sim_idx))

        df = pd.DataFrame(all_data)
        os.makedirs(os.path.dirname(os.path.abspath(file_path)), exist_ok=True)
        df.to_parquet(file_path, index=False)
        logger.info(f"Successfully exported {len(all_data)} rows to {file_path}")

    def export_results_json(self, results: Dict[str, Any], file_path: str):
        """
        Exports results to a JSON format for broader compatibility (SED-ML like structure).
        """
        export_data: Dict[str, Any] = {
            "metadata": {
                "model_name": self.model_name,
                "time_unit": self.time_unit,
                "concentration_unit": self.concentration_unit,
                "basal_growth": results.get("basal_growth"),
                "drug_concentration": self.drug_concentration
            },
            "simulations": []
        }

        for i, hist in enumerate(results.get("histories", [])):
            sim_entry = {
                "id": i,
                "drug_kd": hist["drug_kd"],
                "params": hist.get("params", {}),
                "cell_death": hist.get("cell_death", False),
                "death_step": hist.get("death_step"),
                "data": {
                    "time": hist["time"].tolist(),
                    "growth": hist["growth"].tolist(),
                    "signaling": hist["signaling"].tolist()
                }
            }
            export_data["simulations"].append(sim_entry)

        if "sensitivity" in results:
            export_data["sensitivity"] = results["sensitivity"]

        os.makedirs(os.path.dirname(os.path.abspath(file_path)), exist_ok=True)
        with open(file_path, "w") as f:
            json.dump(export_data, f, indent=2)
        logger.info(f"Successfully exported results to {file_path}")

