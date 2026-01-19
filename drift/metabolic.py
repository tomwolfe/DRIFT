from cobra.io import load_model
import numpy as np
import logging
from typing import List, Dict, Any, Union, Optional, Callable

# Configure logging
logger = logging.getLogger(__name__)


def sigmoidal_mapping(x, k=10, x0=0.5):
    """Sigmoidal mapping function: 1 / (1 + exp(-k * (x - x0)))"""
    return 1 / (1 + np.exp(-k * (x - x0)))


def michaelis_menten_mapping(x, km=0.5):
    """Michaelis-Menten mapping function: x / (km + x)"""
    # Normalize so that at x=1.0 it returns something close to 1.0 (or at least defined)
    # MM is x/(Km+x). At x=1, it is 1/(Km+1).
    # To have it range from 0 to 1, we can use x*(Km+1)/(Km+x)
    return (x * (km + 1)) / (km + x)


class MetabolicBridge:
    """Maps signaling protein concentrations to metabolic Vmax constraints."""

    def __init__(
        self, 
        mappings: Optional[List[Dict[str, Any]]] = None, 
        reverse_mappings: Optional[List[Dict[str, Any]]] = None,
        species_names: Optional[List[str]] = None
    ):
        """
        Initialize the MetabolicBridge.

        Args:
            mappings: List of dicts specifying how signaling affects metabolism.
            reverse_mappings: List of dicts specifying how metabolism affects signaling.
            species_names: Optional list of signaling species names to resolve name-based mappings.
        """
        if mappings is None:
            # Default: mTOR increases glucose uptake capacity
            self.mappings: List[Dict[str, Any]] = [
                {
                    "protein_name": "mTOR",
                    "protein_idx": 2,  # Keep for backward compatibility with default 3-species model
                    "reaction_id": "EX_glc__D_e",
                    "influence": "positive",
                    "base_vmax": 10.0,
                }
            ]
        else:
            self.mappings = mappings
        
        self.reverse_mappings = reverse_mappings or []
        self.species_names = species_names

    def get_feedback(self, fluxes: Dict[str, float]) -> np.ndarray:
        """
        Calculates metabolite-specific feedback signals.
        
        Returns:
            np.ndarray: Vector of feedback signals, one for each signaling species.
        """
        if self.species_names is None:
            # Fallback to global growth-based scalar if no species names defined
            return np.array([self._get_global_feedback(fluxes)])

        feedback_vec = np.ones(len(self.species_names))
        
        # Determine global/energy state
        global_fb = self._get_global_feedback(fluxes)
        energy_fb = self._get_energy_feedback(fluxes)

        # Apply specific reverse mappings
        for rev_map in self.reverse_mappings:
            flux_id = rev_map.get("flux_id")
            species_name = rev_map.get("species_name")
            influence = rev_map.get("influence", "positive")
            weight = rev_map.get("weight", 1.0)
            mapping_type = rev_map.get("mapping_type", "linear")
            mapping_params = rev_map.get("mapping_params", {})

            if flux_id in fluxes and species_name in self.species_names:
                idx = self.species_names.index(species_name)
                # Normalize flux relative to a baseline or max expected (default 10.0)
                baseline = rev_map.get("baseline", 10.0)
                flux_val = np.clip(fluxes[flux_id] / baseline, 0.0, 1.0) 
                
                # Apply mapping function
                if mapping_type == "sigmoidal":
                    k = mapping_params.get("k", 10.0)
                    x0 = mapping_params.get("x0", 0.5)
                    transformed_val = sigmoidal_mapping(flux_val, k=k, x0=x0)
                elif mapping_type == "michaelis-menten":
                    km = mapping_params.get("km", 0.5)
                    transformed_val = michaelis_menten_mapping(flux_val, km=km)
                else: # Default linear
                    transformed_val = flux_val

                if influence == "positive":
                    feedback_vec[idx] = (1.0 - weight) + weight * transformed_val
                else:
                    feedback_vec[idx] = (1.0 - weight) + weight * (1.0 - transformed_val)

        # Apply default biological feedbacks if not explicitly overridden
        for i, name in enumerate(self.species_names):
            if feedback_vec[i] == 1.0: # Only if not already set by specific mapping
                if name == "mTOR":
                    # mTOR is sensitive to both growth and energy
                    feedback_vec[i] = 0.5 * global_fb + 0.5 * energy_fb
                elif name == "AMPK":
                    # AMPK is activated (high value) when energy is low (low energy_fb)
                    # So AMPK feedback = 1.0 - energy_fb
                    feedback_vec[i] = 1.0 - energy_fb
                elif name in ["PI3K", "AKT"]:
                    # These might have some mild sensitivity to global state
                    feedback_vec[i] = 0.8 + 0.2 * global_fb

        return feedback_vec

    def _get_global_feedback(self, fluxes: Dict[str, float]) -> float:
        """Calculates a global feedback signal based on growth rate (internal helper)."""
        if not fluxes:
            return 1.0
            
        growth_keys = [
            "Biomass_Ecoli_core", 
            "BIOMASS_Ecoli_core_w_GAM", 
            "growth", 
            "BIOMASS_RECON1", 
            "BIOMASS_reaction",
            "BIOMASS_maintenance"
        ]
        growth_flux = 0.0
        for key in growth_keys:
            if key in fluxes:
                growth_flux = fluxes[key]
                break
        
        if growth_flux == 0.0 and fluxes:
            # Fallback to any positive biomass-like flux if common names not found
            growth_flux = max(fluxes.get(k, 0.0) for k in fluxes if "BIOMASS" in k.upper()) if any("BIOMASS" in k.upper() for k in fluxes) else 0.0

        feedback = float(np.clip(growth_flux / 0.2, 0.0, 1.0))
        return feedback

    def _get_energy_feedback(self, fluxes: Dict[str, float]) -> float:
        """Calculates an energy-sensing feedback signal (ATP/ADP ratio proxy)."""
        if not fluxes:
            return 1.0
            
        # Proxy: Total ATP production flux
        # In textbook model, ATPS4r is the main ATP synthase
        atp_flux = fluxes.get("ATPS4r", 0.0)
        
        # If not found, look for other ATP-producing reactions
        if atp_flux == 0.0:
            atp_flux = sum(f for k, f in fluxes.items() if "ATP" in k and f > 0)
            
        # Normalize: textbook ATPS4r is around 2-10 under good growth
        energy_state = float(np.clip(atp_flux / 5.0, 0.0, 1.0))
        return energy_state

    def get_constraints(self, signaling_state: Union[List[float], np.ndarray]) -> Dict[str, float]:
        """
        Translates signaling state to a dict of {reaction_id: lower_bound}.

        Args:
            signaling_state: Normalized protein concentrations

        Returns:
            dict: Dictionary mapping reaction IDs to constraint values
        """
        if (
            not isinstance(signaling_state, (list, tuple, np.ndarray))
        ):
            raise ValueError(
                f"signaling_state must be a list, tuple, or array, got {type(signaling_state)}"
            )

        constraints: Dict[str, float] = {}
        state_len = len(signaling_state)

        for map_config in self.mappings:
            idx = map_config.get("protein_idx")
            name = map_config.get("protein_name")
            rxn_id: str = map_config["reaction_id"]
            base_vmax: float = map_config.get("base_vmax", 10.0)
            influence: str = map_config.get("influence", "positive")
            mapping_fn = map_config.get("mapping_fn")

            # Resolve name to index if possible
            if name is not None and self.species_names is not None:
                if name in self.species_names:
                    idx = self.species_names.index(name)
                else:
                    logger.warning(f"Protein name '{name}' not found in species_names {self.species_names}")

            if idx is None or not (0 <= idx < state_len):
                logger.warning(f"Invalid protein index {idx} (name: {name}) for state of length {state_len}, skipping mapping")
                continue

            protein_level = float(signaling_state[idx])

            if not (0 <= protein_level <= 1):
                logger.warning(
                    f"Protein level {protein_level} out of range [0,1], clamping"
                )
                protein_level = max(0.0, min(1.0, protein_level))

            if mapping_fn is not None:
                # Use custom mapping function
                scaling = mapping_fn(protein_level)
            elif influence == "positive":
                # Default linear scaling: 0.1 basal to 1.0 max
                scaling = 0.1 + 0.9 * protein_level
            else:
                # Default inhibitory effect
                scaling = 1.0 - 0.9 * protein_level

            # COBRA exchange reactions: Flux > -UB (uptake)
            # We scale the magnitude of the negative lower bound.
            constraints[rxn_id] = -(base_vmax * scaling)

        return constraints


    @classmethod
    def get_human_cancer_bridge(cls):
        """
        Returns a pre-configured MetabolicBridge for human cancer research.
        Maps PI3K/AKT/mTOR signaling to Recon1 metabolic subsystems.
        """
        builder = BridgeBuilder()
        builder.set_species_names(["PI3K", "AKT", "mTOR", "AMPK"])
        
        # Warburg Effect: mTOR increases glucose uptake and lactate secretion
        builder.add_mapping(protein_name="mTOR", reaction_id="EX_glc__D_e", influence="positive", base_vmax=15.0)
        builder.add_mapping(protein_name="mTOR", reaction_id="EX_lac__L_e", influence="positive", base_vmax=20.0)
        
        # Glutaminolysis: AKT increases glutamine uptake
        builder.add_mapping(protein_name="AKT", reaction_id="EX_gln__L_e", influence="positive", base_vmax=5.0)
        
        # Energy sensing: AMPK activates fatty acid oxidation and oxidative phosphorylation
        builder.add_mapping(protein_name="AMPK", reaction_id="EX_o2_e", influence="positive", base_vmax=20.0)
        
        # Reverse mappings (Metabolism -> Signaling)
        # Low energy (low ATP synthase flux) activates AMPK
        builder.add_reverse_mapping(flux_id="ATPS4r", species_name="AMPK", influence="negative", weight=0.8)
        
        # Growth influences mTOR
        builder.add_reverse_mapping(flux_id="BIOMASS_RECON1", species_name="mTOR", influence="positive", weight=0.5)
        
        return builder.build()


class BridgeBuilder:
    """Fluent API for building MetabolicBridge instances."""

    def __init__(self):
        self.mappings = []
        self.reverse_mappings = []
        self.species_names = None

    def set_species_names(self, species_names: List[str]) -> "BridgeBuilder":
        """Sets the species names for the bridge."""
        self.species_names = species_names
        return self

    def add_mapping(
        self, 
        protein: Optional[Union[int, str]] = None, 
        reaction_id: str = "", 
        influence: str = "positive", 
        base_vmax: float = 10.0,
        mapping_fn: Optional[Callable[[float], float]] = None,
        protein_idx: Optional[int] = None,
        protein_name: Optional[str] = None
    ) -> "BridgeBuilder":
        """
        Adds a mapping from a signaling protein to a metabolic reaction.
        """
        if influence not in ["positive", "negative"]:
            raise ValueError("influence must be 'positive' or 'negative'")
        
        mapping = {
            "reaction_id": reaction_id,
            "influence": influence,
            "base_vmax": base_vmax,
            "mapping_fn": mapping_fn
        }
        
        # Resolve protein identifier
        p = protein
        if p is None:
            p = protein_idx
        if p is None:
            p = protein_name
            
        if p is None:
            raise ValueError("Must provide protein index or name")
        
        if isinstance(p, int):
            mapping["protein_idx"] = p
        else:
            mapping["protein_name"] = p

        self.mappings.append(mapping)
        return self

    def add_reverse_mapping(
        self,
        flux_id: str,
        species_name: str,
        influence: str = "positive",
        weight: float = 1.0,
        mapping_type: str = "linear",
        mapping_params: Optional[Dict[str, Any]] = None,
        baseline: float = 10.0
    ) -> "BridgeBuilder":
        """
        Adds a feedback mapping from a metabolic flux to a signaling species.
        """
        self.reverse_mappings.append({
            "flux_id": flux_id,
            "species_name": species_name,
            "influence": influence,
            "weight": weight,
            "mapping_type": mapping_type,
            "mapping_params": mapping_params or {},
            "baseline": baseline
        })
        return self

    def build(self) -> MetabolicBridge:
        """Returns the constructed MetabolicBridge."""
        return MetabolicBridge(
            mappings=self.mappings, 
            reverse_mappings=self.reverse_mappings,
            species_names=self.species_names
        )


class DFBASolver:
    """Dynamic Flux Balance Analysis solver."""

    def __init__(self, model_name="textbook"):
        """
        Initialize the DFBASolver.

        Args:
            model_name (str): Name of the metabolic model to use
        """
        self.model_name = model_name
        self.available_solvers = self._check_solver()
        
        # Headless mode: If no solvers, don't attempt to load or validate
        if not self.available_solvers:
            logger.warning(f"DFBASolver initialized in HEADLESS mode (no solver found). FBA will be disabled.")
            self.model = None
            return

        self.model = self._load_model_safe(model_name)
        
        # Pareto: Immediate validation to ensure the model is actually usable.
        if self.model:
            self.validate_model()

    def validate_model(self):
        """
        Performs a pre-flight check to ensure the model is viable for DRIFT.
        Checks for solvers, objectives, and basic growth capability.
        """
        if not self.model:
            return

        logger.info(f"[*] Validating model '{self.model_name}' for DRIFT compatibility...")
        
        # 1. Check for objective
        if not self.model.objective:
            raise RuntimeError(f"Model '{self.model_name}' has no defined objective function.")
            
        # 2. Check if model is solvable in default state
        try:
            solution = self.model.optimize()
            if solution.status != "optimal":
                raise RuntimeError(f"Model '{self.model_name}' is infeasible in its default state (Status: {solution.status}).")
            if solution.objective_value <= 0:
                logger.warning(f"Model '{self.model_name}' has zero growth in default state. Check medium/bounds.")
        except Exception as e:
            raise RuntimeError(f"Critical error during model validation: {str(e)}")

        logger.info(f"[+] Model '{self.model_name}' validated successfully.")

    def _check_solver(self):
        """Checks if a valid COBRA solver is available and returns a list of them."""
        try:
            from optlang import available_solvers
            active_solvers = [s for s, available in available_solvers.items() if available]
            
            if not active_solvers:
                error_msg = (
                    "No COBRA-compatible solvers (GLPK, CPLEX, GUROBI, etc.) found.\n"
                    "FBA requires a linear programming solver to function.\n"
                    "Recommended: Install GLPK by running 'pip install swiglpk' or 'conda install glpk'."
                )
                logger.warning(error_msg)
                # We no longer raise here to allow "headless" operation for signaling-only runs
                return []
                
            logger.info(f"Available solvers: {active_solvers}. Using default.")
            return active_solvers
        except ImportError:
            logger.warning("Could not check for available solvers via optlang. FBA might fail.")
            return []

    def check_solver_sensitivity(self):
        """
        Runs the current model with all available solvers and compares results.
        This is critical for ensuring reproducibility across different environments.
        """
        if not self.model or not self.available_solvers:
            logger.warning("Solver sensitivity check skipped: No model or solver available.")
            return {}

        import time
        from optlang import available_solvers
        
        active_solvers = [s for s, available in available_solvers.items() if available]
        logger.info(f"[*] Running solver sensitivity check across: {active_solvers}")
        
        results = {}
        original_solver = self.model.solver.interface.__name__.split('.')[-1].replace('interface', '').lower()
        
        for solver_name in active_solvers:
            try:
                # Try setting the solver (cobra often accepts lowercase)
                try:
                    self.model.solver = solver_name
                except Exception:
                    try:
                        self.model.solver = solver_name.lower()
                    except Exception:
                        logger.debug(f"Could not set solver to {solver_name}, skipping.")
                        continue
                    
                start_time = time.time()
                solution = self.model.optimize()
                elapsed = time.time() - start_time
                
                results[solver_name] = {
                    "objective": solution.objective_value,
                    "status": solution.status,
                    "time": elapsed
                }
                logger.info(f"  - {solver_name}: obj={solution.objective_value:.6f}, time={elapsed:.4f}s")
            except Exception as e:
                logger.warning(f"  - {solver_name}: FAILED ({str(e)})")
        
        # Restore original solver
        try:
            self.model.solver = original_solver
        except Exception as e:
            logger.debug(f"Could not restore original solver {original_solver}: {e}")

        # Check for discrepancies
        if len(results) > 1:
            objectives = [r["objective"] for r in results.values() if r["status"] == "optimal"]
            if objectives:
                max_diff = max(objectives) - min(objectives)
                if max_diff > 1e-6:
                    logger.warning(f"[!] Significant solver sensitivity detected! Max difference in objective: {max_diff:.8f}")
                else:
                    logger.info("[+] Solver results are consistent (diff < 1e-6).")
        
        return results

    def _load_model_safe(self, name):
        """
        Safely load a metabolic model with fallback options.
        """
        if not self.available_solvers:
            return None

        try:
            model = load_model(name)
            logger.info(f"Successfully loaded model: {name}")
            return model
        except Exception as e:
            print(f"[!] Warning: Failed to load requested model '{name}'.")
            logger.warning(
                f"Failed to load model '{name}': {str(e)}. Attempting fallbacks..."
            )

            fallback_models = ["textbook", "e_coli_core", "iJO1366", "recon1"]
            if name in fallback_models:
                fallback_models.remove(name)

            for fallback_name in fallback_models:
                try:
                    logger.info(f"Trying fallback model: {fallback_name}")
                    model = load_model(fallback_name)
                    print(f"[!] Falling back to model: '{fallback_name}'. Results may differ from expectations.")
                    logger.info(f"Successfully loaded fallback model: {fallback_name}")
                    return model
                except Exception as fallback_error:
                    logger.warning(
                        f"Failed to load fallback model '{fallback_name}': {str(fallback_error)}"
                    )
                    continue

            error_msg = (
                f"Could not load requested model '{name}' or any fallback models. "
                "Please ensure the model name is valid and the model is available."
            )
            logger.critical(error_msg)
            raise RuntimeError(error_msg) from e

    def solve_step(self, constraints):
        """
        Solves FBA with specific constraints and provides diagnostics on failure.
        """
        if not self.model:
            # Headless mode: return a dummy result
            return {
                "objective_value": 0.0,
                "fluxes": {},
                "status": "headless",
                "diagnostic": "No solver available (Headless Mode)"
            }

        if not isinstance(constraints, dict):
            raise ValueError(f"constraints must be a dictionary, got {type(constraints)}")

        # Track applied constraints for diagnostics
        applied_constraints = {}
        for rxn_id, lb in constraints.items():
            try:
                rxn = self.model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = lb
                applied_constraints[rxn_id] = lb
            except KeyError:
                logger.debug(f"Reaction {rxn_id} not found in model {self.model_name}")
                continue

        try:
            solution = self.model.optimize()
            
            if solution.status != "optimal":
                # Pareto: Diagnostic "Autopsy"
                # Identify limiting factors by checking which applied constraints are at their bounds
                limiting_factors = []
                for rxn_id, lb in applied_constraints.items():
                    if lb > -1.0: # If we've significantly restricted an uptake
                        limiting_factors.append(f"{rxn_id} (LB: {lb:.2f})")
                
                # If still unknown, look for the 'bottleneck' via relaxation (Pareto approach)
                diag_msg = "Global Infeasibility"
                if limiting_factors:
                    diag_msg = f"Potential bottlenecks: {', '.join(limiting_factors[:3])}"
                
                # Try to find the specific metabolite that is limiting
                # In FBA, this is often the one where the shadow price would be highest if it were feasible
                # But since it's infeasible, we can look at the 'irreducible inconsistent subsystem' (IIS)
                # For now, we'll provide the list of restrictive constraints.
                
                return {
                    "objective_value": 0.0,
                    "fluxes": {},
                    "status": solution.status,
                    "diagnostic": diag_msg,
                    "bottlenecks": applied_constraints
                }

            return {
                "objective_value": solution.objective_value,
                "fluxes": solution.fluxes.to_dict(),
                "status": solution.status
            }
        except Exception as e:
            logger.error(f"FBA optimization failed: {e}")
            return {
                "objective_value": 0.0,
                "fluxes": {},
                "status": "failed",
                "error": str(e),
                "diagnostic": "Solver error or numerical instability"
            }