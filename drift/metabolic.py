from cobra.io import load_model
import numpy as np
import logging
from typing import List, Dict, Any, Union, Optional

# Configure logging
logger = logging.getLogger(__name__)


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

    def get_feedback(self, fluxes: Dict[str, float]) -> float:
        """
        Calculates a global feedback signal based on metabolic state.
        By default, it uses the growth rate (objective) normalized to a typical value.
        """
        if not fluxes:
            return 1.0  # Default to full activity if no fluxes provided
            
        # Try to find a growth-related flux
        growth_keys = ["Biomass_Ecoli_core", "BIOMASS_Ecoli_core_w_GAM", "growth"]
        growth_flux = 0.0
        for key in growth_keys:
            if key in fluxes:
                growth_flux = fluxes[key]
                break
        
        if growth_flux == 0.0 and fluxes:
            # If no explicit growth key, use the max flux as a proxy for activity
            growth_flux = max(fluxes.values()) if fluxes else 0.0

        # Normalize: 0.2 is a typical max growth for core models
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

            if influence == "positive":
                # Scaling factor: 0.1 basal to 1.0 max
                scaling = 0.1 + 0.9 * protein_level
            else:
                # Inhibitory effect
                scaling = 1.0 - 0.9 * protein_level

            # COBRA exchange reactions: Flux > -UB (uptake)
            # We scale the magnitude of the negative lower bound.
            constraints[rxn_id] = -(base_vmax * scaling)

        return constraints


class BridgeBuilder:
    """Fluent API for building MetabolicBridge instances."""

    def __init__(self):
        self.mappings = []
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
        protein_idx: Optional[int] = None,
        protein_name: Optional[str] = None
    ) -> "BridgeBuilder":
        """
        Adds a mapping from a signaling protein to a metabolic reaction.
        
        Args:
            protein: Either the index (int) or the name (str) of the signaling protein.
            reaction_id: The ID of the metabolic reaction to constrain.
            influence: 'positive' or 'negative' effect.
            base_vmax: The maximum flux capacity.
            protein_idx: Explicit protein index (for backward compatibility).
            protein_name: Explicit protein name.
        """
        if influence not in ["positive", "negative"]:
            raise ValueError("influence must be 'positive' or 'negative'")
        
        mapping = {
            "reaction_id": reaction_id,
            "influence": influence,
            "base_vmax": base_vmax
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

    def build(self) -> MetabolicBridge:
        """Returns the constructed MetabolicBridge."""
        return MetabolicBridge(mappings=self.mappings, species_names=self.species_names)


class DFBASolver:
    """Dynamic Flux Balance Analysis solver."""

    def __init__(self, model_name="textbook"):
        """
        Initialize the DFBASolver.

        Args:
            model_name (str): Name of the metabolic model to use
        """
        self.model_name = model_name
        self.model = self._load_model_safe(model_name)

    def _load_model_safe(self, name):
        """
        Safely load a metabolic model with fallback options.

        Args:
            name (str): Name of the model to load

        Returns:
            cobra.Model: Loaded metabolic model
        """
        # First try the requested model
        try:
            model = load_model(name)
            logger.info(f"Successfully loaded model: {name}")
            return model
        except Exception as e:
            logger.warning(
                f"Failed to load model '{name}': {str(e)}. Attempting fallbacks..."
            )

            # Try common model names as fallbacks
            fallback_models = ["e_coli_core", "textbook"]
            if name in fallback_models:
                fallback_models.remove(name)  # Don't retry the same failed name

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

            # If all attempts fail, raise a more informative error
            error_msg = (
                f"Could not load requested model '{name}' or any fallback models. "
                f"Please ensure the model name is valid and the model is available. "
                f"Common valid models include: 'textbook', 'e_coli_core'."
            )
            logger.critical(error_msg)
            raise RuntimeError(error_msg) from e

    def solve_step(self, constraints):
        """
        Solves FBA with specific constraints.

        Args:
            constraints: {rxn_id: lower_bound_value}

        Returns:
            tuple: (objective_value, flux_dictionary)
        """
        if not isinstance(constraints, dict):
            raise ValueError(
                f"constraints must be a dictionary, got {type(constraints)}"
            )

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
