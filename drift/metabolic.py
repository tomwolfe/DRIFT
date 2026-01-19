import cobra
from cobra.io import load_model
import numpy as np
import logging

# Configure logging
logger = logging.getLogger(__name__)

class MetabolicBridge:
    """Maps signaling protein concentrations to metabolic Vmax constraints."""
    def __init__(self, mappings=None):
        """
        mappings: List of dicts specifying how signaling affects metabolism.
        Example: [
            {'protein_idx': 2, 'reaction_id': 'EX_glc__D_e', 'influence': 'positive', 'base_vmax': 10.0}
        ]
        """
        if mappings is None:
            # Default: mTOR (idx 2) increases glucose uptake capacity
            self.mappings = [
                {'protein_idx': 2, 'reaction_id': 'EX_glc__D_e', 'influence': 'positive', 'base_vmax': 10.0}
            ]
        else:
            self.mappings = mappings

    def get_constraints(self, signaling_state):
        """
        Translates signaling state to a dict of {reaction_id: lower_bound}.
        signaling_state: [PI3K, AKT, mTOR]
        """
        constraints = {}
        
        for map_config in self.mappings:
            idx = map_config['protein_idx']
            rxn_id = map_config['reaction_id']
            base_vmax = map_config.get('base_vmax', 10.0)
            influence = map_config.get('influence', 'positive')
            
            protein_level = signaling_state[idx]
            
            if influence == 'positive':
                # Scaling factor: 0.1 basal to 1.0 max
                scaling = 0.1 + 0.9 * protein_level
            else:
                # Inhibitory effect
                scaling = 1.0 - 0.9 * protein_level
                
            # COBRA exchange reactions: Flux > -UB (uptake)
            # We scale the magnitude of the negative lower bound.
            constraints[rxn_id] = - (base_vmax * scaling)
            
        return constraints

class DFBASolver:
    def __init__(self, model_name='textbook'):
        self.model_name = model_name
        self.model = self._load_model_safe(model_name)
        
    def _load_model_safe(self, name):
        try:
            return load_model(name)
        except Exception as e:
            logger.warning(f"Failed to load model '{name}': {e}. Falling back to 'e_coli_core'.")
            try:
                return load_model('e_coli_core')
            except Exception as e2:
                logger.error(f"Critical error: Could not load fallback model: {e2}")
                raise

    def solve_step(self, constraints):
        """
        Solves FBA with specific constraints.
        constraints: {rxn_id: lower_bound_value}
        """
        for rxn_id, lb in constraints.items():
            try:
                rxn = self.model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = lb
            except KeyError:
                logger.debug(f"Reaction {rxn_id} not found in model {self.model_name}")
                continue
        
        try:
            solution = self.model.optimize()
            if solution.status == 'optimal':
                return solution.objective_value, solution.fluxes.to_dict()
        except Exception as e:
            logger.error(f"FBA optimization failed: {e}")
            
        return 0.0, {}
