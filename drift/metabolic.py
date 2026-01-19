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


def fuzzy_match_reaction_id(query_id: str, model_reaction_ids: List[str]) -> Optional[str]:
    """
    Attempts to find a matching reaction ID in the model using common variations.
    Handles 'EX_' prefixes, case sensitivity, and double underscores.
    """
    if query_id in model_reaction_ids:
        return query_id
    
    # Try case-insensitive
    lower_ids = {rid.lower(): rid for rid in model_reaction_ids}
    if query_id.lower() in lower_ids:
        return lower_ids[query_id.lower()]
    
    # Try adding/removing 'EX_' prefix
    variations = []
    if query_id.startswith("EX_"):
        variations.append(query_id[3:])
    else:
        variations.append(f"EX_{query_id}")
    
    # Try replacing double underscores with single or vice-versa
    if "__" in query_id:
        variations.append(query_id.replace("__", "_"))
    else:
        # This is risky but sometimes helpful
        pass

    for var in variations:
        if var in model_reaction_ids:
            return var
        if var.lower() in lower_ids:
            return lower_ids[var.lower()]
            
    return None


class MetabolicBridge:
    """Maps signaling protein concentrations to metabolic Vmax constraints."""

    def __init__(
        self, 
        mappings: Optional[List[Dict[str, Any]]] = None, 
        reverse_mappings: Optional[List[Dict[str, Any]]] = None,
        species_names: Optional[List[str]] = None,
        strict_mapping: bool = False
    ):
        """
        Initialize the MetabolicBridge.

        Args:
            mappings: List of dicts specifying how signaling affects metabolism.
            reverse_mappings: List of dicts specifying how metabolism affects signaling.
            species_names: Optional list of signaling species names to resolve name-based mappings.
            strict_mapping: If True, disables fuzzy matching and raises error on missing IDs.
        """
        self.strict_mapping = strict_mapping
        if mappings is None:
            # Default: mTOR increases glucose uptake capacity
            self.mappings: List[Dict[str, Any]] = [
                {
                    "protein_name": "mTOR",
                    "reaction_id": "EX_glc__D_e",
                    "influence": "positive",
                    "base_vmax": 10.0,
                    "basal_scaling": 0.1,
                    "max_scaling": 0.9
                }
            ]
            # Provide default species names if they are missing to resolve "mTOR"
            if species_names is None:
                self.species_names: Optional[List[str]] = ["PI3K", "AKT", "mTOR"]
            else:
                self.species_names = species_names
        else:
            self.mappings = mappings
            self.species_names = species_names
        
        self.reverse_mappings = reverse_mappings or []

    def validate_with_model(self, model: Any) -> bool:
        """
        Validates and repairs reaction IDs in mappings using the provided model.
        Uses fuzzy matching to resolve common naming discrepancies unless strict_mapping is True.
        
        Returns:
            bool: True if all mappings were resolved, False otherwise.
        """
        if model is None:
            return False
            
        model_rxn_ids = [r.id for r in model.reactions]
        all_resolved = True
        
        # Validate forward mappings
        for mapping in self.mappings:
            rxn_id = mapping.get("reaction_id")
            if rxn_id not in model_rxn_ids:
                if self.strict_mapping:
                    logger.error(f"STRICT MAPPING FAILURE: Reaction ID '{rxn_id}' not found in model.")
                    all_resolved = False
                    continue

                matched = fuzzy_match_reaction_id(rxn_id, model_rxn_ids)
                if matched:
                    logger.warning(f"CRITICAL: Fuzzy matched mapping reaction '{rxn_id}' to '{matched}'. "
                                 "This may lead to biological misassignment. Use strict_mapping=True to prevent this.")
                    mapping["reaction_id"] = matched
                else:
                    logger.warning(f"Could not resolve reaction ID '{rxn_id}' in model.")
                    all_resolved = False
                    
        # Validate reverse mappings
        for rev_map in self.reverse_mappings:
            flux_id = rev_map.get("flux_id")
            if flux_id not in model_rxn_ids:
                if self.strict_mapping:
                    logger.error(f"STRICT MAPPING FAILURE: Reverse flux ID '{flux_id}' not found in model.")
                    all_resolved = False
                    continue

                matched = fuzzy_match_reaction_id(flux_id, model_rxn_ids)
                if matched:
                    logger.warning(f"CRITICAL: Fuzzy matched reverse mapping flux '{flux_id}' to '{matched}'.")
                    rev_map["flux_id"] = matched
                else:
                    logger.warning(f"Could not resolve reverse mapping flux ID '{flux_id}' in model.")
                    all_resolved = False
                    
        return all_resolved

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

        # If no specific mappings were applied, apply a generic global feedback
        # to ensure the system remains coupled even with minimal configuration.
        for i in range(len(feedback_vec)):
            if feedback_vec[i] == 1.0:
                # Default behavior: slight sensitivity to global metabolic state
                feedback_vec[i] = 0.9 + 0.1 * global_fb

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
        # ATPS4r (textbook), ATPMS (Recon1), etc.
        atp_keys = ["ATPS4r", "ATPMS", "ATPS", "ATPSv"]
        atp_flux = 0.0
        for key in atp_keys:
            if key in fluxes:
                atp_flux = fluxes[key]
                break
        
        # If not found, look for other ATP-producing reactions
        if atp_flux == 0.0:
            # Sum of all positive fluxes through reactions containing 'ATP' and 'synthase' or 'synth'
            atp_flux = sum(f for k, f in fluxes.items() if ("ATP" in k.upper() and ("SYNTH" in k.upper() or "S4R" in k.upper())) and f > 0)
            
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

            # Resolve name to index (Enforced Default)
            if name is not None:
                if self.species_names is not None:
                    if name in self.species_names:
                        idx = self.species_names.index(name)
                    else:
                        raise ValueError(f"Protein name '{name}' not found in species_names {self.species_names}. "
                                         "Ensure your Topology matches your MetabolicBridge.")
                elif idx is None:
                    # No species names and no index provided
                    raise ValueError(f"Protein name '{name}' provided but MetabolicBridge has no 'species_names' "
                                     "to resolve it, and no 'protein_idx' was specified.")

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
            else:
                # Configuration-based scaling (Pareto improvement: avoid hardcoded arbitrary constants)
                basal = map_config.get("basal_scaling", 0.1)
                max_s = map_config.get("max_scaling", 0.9)
                
                if influence == "positive":
                    scaling = basal + max_s * protein_level
                else:
                    scaling = (basal + max_s) - max_s * protein_level

            # COBRA exchange reactions: Flux > -UB (uptake)
            # scale the magnitude of the negative lower bound.
            constraints[rxn_id] = -(base_vmax * scaling)

        return constraints


class BridgeBuilder:
    """Fluent API for building MetabolicBridge instances."""

    def __init__(self):
        self.mappings = []
        self.reverse_mappings = []
        self.species_names: Optional[List[str]] = None
        self.strict_mapping = False

    def set_strict_mapping(self, strict: bool = True) -> "BridgeBuilder":
        """Enables or disables strict reaction ID mapping."""
        self.strict_mapping = strict
        return self

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
        protein_name: Optional[str] = None,
        basal_scaling: float = 0.1,
        max_scaling: float = 0.9
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
            "mapping_fn": mapping_fn,
            "basal_scaling": basal_scaling,
            "max_scaling": max_scaling
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
            species_names=self.species_names,
            strict_mapping=self.strict_mapping
        )


class DFBASolver:
    """Dynamic Flux Balance Analysis solver."""

    def __init__(self, model_name="textbook", strict=False):
        """
        Initialize the DFBASolver.

        Args:
            model_name (str): Name of the metabolic model to use
            strict (bool): If True, raises an error if no LP solver is found.
        """
        self.model_name = model_name
        self.strict = strict
        self.available_solvers = self._check_solver()
        self.headless = not self.available_solvers
        
        # Headless mode: If no solvers, don't attempt to load or validate
        if self.headless:
            if self.strict:
                raise RuntimeError(
                    "DFBASolver: No LP solver found and 'strict' mode is enabled. "
                    "Install GLPK (pip install swiglpk) or another COBRA-compatible solver."
                )
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
            if self.strict:
                raise RuntimeError("DFBASolver: Attempted to solve FBA in strict mode with no available solver.")
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