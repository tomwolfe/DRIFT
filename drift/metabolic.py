from cobra.io import load_model
import numpy as np
import logging
from typing import List, Dict, Any, Union, Optional

# Configure logging
logger = logging.getLogger(__name__)


class MetabolicBridge:
    """Maps signaling protein concentrations to metabolic Vmax constraints."""

    def __init__(self, mappings: Optional[List[Dict[str, Any]]] = None):
        """
        Initialize the MetabolicBridge.

        Args:
            mappings: List of dicts specifying how signaling affects metabolism.
                     Example: [
                         {'protein_idx': 2, 'reaction_id': 'EX_glc__D_e', 'influence': 'positive', 'base_vmax': 10.0}
                     ]
        """
        if mappings is None:
            # Default: mTOR (idx 2) increases glucose uptake capacity
            self.mappings: List[Dict[str, Any]] = [
                {
                    "protein_idx": 2,
                    "reaction_id": "EX_glc__D_e",
                    "influence": "positive",
                    "base_vmax": 10.0,
                }
            ]
        else:
            self.mappings = mappings

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
            idx: int = map_config["protein_idx"]
            rxn_id: str = map_config["reaction_id"]
            base_vmax: float = map_config.get("base_vmax", 10.0)
            influence: str = map_config.get("influence", "positive")

            if not (0 <= idx < state_len):
                logger.warning(f"Invalid protein index {idx} for state of length {state_len}, skipping mapping")
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
            if solution.status == "optimal":
                return solution.objective_value, solution.fluxes.to_dict()
            else:
                logger.warning(
                    f"FBA optimization returned non-optimal status: {solution.status}"
                )
                return 0.0, {}
        except Exception as e:
            logger.error(f"FBA optimization failed: {e}")
            return 0.0, {}
