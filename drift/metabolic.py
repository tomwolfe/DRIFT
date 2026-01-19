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
                    "protein_idx": 2,  # Keep for backward compatibility
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
        global_fb = self._get_global_feedback(fluxes)

        # Apply specific reverse mappings
        for rev_map in self.reverse_mappings:
            flux_id = rev_map.get("flux_id")
            species_name = rev_map.get("species_name")
            influence = rev_map.get("influence", "positive")
            weight = rev_map.get("weight", 1.0)

            if flux_id in fluxes and species_name in self.species_names:
                idx = self.species_names.index(species_name)
                # Normalize flux (assuming 0 to 1 range for simplified models, or relative to a baseline)
                # For this implementation, we'll use a simple linear mapping for demonstration
                flux_val = np.clip(fluxes[flux_id] / 10.0, 0.0, 1.0) 
                
                if influence == "positive":
                    feedback_vec[idx] = (1.0 - weight) + weight * flux_val
                else:
                    feedback_vec[idx] = (1.0 - weight) + weight * (1.0 - flux_val)

        # For mTOR (index 2 in default), if no specific feedback, use global
        if "mTOR" in self.species_names:
            idx = self.species_names.index("mTOR")
            # Only apply global if specific feedback hasn't been set (remains 1.0)
            if feedback_vec[idx] == 1.0:
                feedback_vec[idx] = global_fb
        elif len(feedback_vec) > 2: # Backward compatibility for 3-species model
            if feedback_vec[2] == 1.0:
                feedback_vec[2] = global_fb

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
            growth_flux = max(fluxes.values()) if fluxes else 0.0

        feedback = float(np.clip(growth_flux / 0.2, 0.0, 1.0))
        return feedback

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
        weight: float = 1.0
    ) -> "BridgeBuilder":
        """
        Adds a feedback mapping from a metabolic flux to a signaling species.
        """
        self.reverse_mappings.append({
            "flux_id": flux_id,
            "species_name": species_name,
            "influence": influence,
            "weight": weight
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
        self._check_solver()
        self.model = self._load_model_safe(model_name)

    def _check_solver(self):
        """Checks if a valid COBRA solver is available and raises a clear error if not."""
        try:
            from optlang import available_solvers
            active_solvers = [s for s, available in available_solvers.items() if available]
            
            if not active_solvers:
                error_msg = (
                    "No COBRA-compatible solvers (GLPK, CPLEX, GUROBI, etc.) found.\n"
                    "FBA requires a linear programming solver to function.\n"
                    "Recommended: Install GLPK by running 'pip install swiglpk' or 'conda install glpk'."
                )
                logger.critical(error_msg)
                raise RuntimeError(error_msg)
                
            logger.info(f"Available solvers: {active_solvers}. Using default.")
        except ImportError:
            logger.warning("Could not check for available solvers via optlang. FBA might fail.")

    def _load_model_safe(self, name):
        """
        Safely load a metabolic model with fallback options.
        """
        try:
            model = load_model(name)
            logger.info(f"Successfully loaded model: {name}")
            return model
        except Exception as e:
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
        Solves FBA with specific constraints.
        """
        if not isinstance(constraints, dict):
            raise ValueError(f"constraints must be a dictionary, got {type(constraints)}")

        for rxn_id, lb in constraints.items():
            try:
                rxn = self.model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = lb
            except KeyError:
                logger.debug(f"Reaction {rxn_id} not found in model {self.model_name}")
                continue

        try:
            solution = self.model.optimize()
            return {
                "objective_value": solution.objective_value if solution.status == "optimal" else 0.0,
                "fluxes": solution.fluxes.to_dict() if solution.status == "optimal" else {},
                "status": solution.status
            }
        except Exception as e:
            logger.error(f"FBA optimization failed: {e}")
            return {
                "objective_value": 0.0,
                "fluxes": {},
                "status": "failed",
                "error": str(e)
            }