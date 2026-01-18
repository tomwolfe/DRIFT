import cobra
from cobra.io import load_model
import numpy as np

class MetabolicBridge:
    """Maps signaling protein concentrations to metabolic Vmax constraints."""
    def __init__(self, reaction_map=None):
        # Default map: mTOR activity scales Glucose uptake
        # For 'textbook' model, reaction is 'EX_glc__D_e' (exchange) 
        # but in cobra textbook it's often 'GLCpts' for transport.
        # Let's check common IDs or allow custom mapping.
        if reaction_map is None:
            self.reaction_map = {
                'mtor': 'EX_glc__D_e' 
            }
        else:
            self.reaction_map = reaction_map

    def get_constraints(self, signaling_state, base_vmax=10.0):
        """
        Translates signaling state to a dict of {reaction_id: upper_bound}.
        signaling_state: [PI3K, AKT, mTOR]
        """
        mtor_level = signaling_state[2]
        
        # Heuristic: mTOR increases glucose uptake capacity
        # We model this by adjusting the lower bound of the exchange (for uptake)
        # In COBRA, exchange reactions are usually: Flux > -UB (uptake)
        # So we scale the magnitude of the negative lower bound.
        
        # Scaling factor: 0.1 basal to 1.0 max
        scaling = 0.1 + 0.9 * mtor_level
        constraints = {}
        
        if 'mtor' in self.reaction_map:
            rxn_id = self.reaction_map['mtor']
            constraints[rxn_id] = - (base_vmax * scaling) # Lower bound for uptake
            
        return constraints

class DFBASolver:
    def __init__(self, model_name='textbook'):
        try:
            self.model = load_model(model_name)
        except:
            # Fallback to e_coli_core if textbook fails
            self.model = load_model('e_coli_core')
            
    def solve_step(self, constraints):
        """
        Solves FBA with specific constraints.
        constraints: {rxn_id: lower_bound_value}
        """
        with self.model as m:
            for rxn_id, lb in constraints.items():
                if rxn_id in m.reactions:
                    # For exchange reactions, uptake is usually the lower bound (negative)
                    m.reactions.get_by_id(rxn_id).lower_bound = lb
            
            solution = m.optimize()
            if solution.status == 'optimal':
                return solution.objective_value, solution.fluxes.to_dict()
            else:
                return 0.0, {}
