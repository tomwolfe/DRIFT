import numpy as np
import json
import logging
from typing import List, Dict, Any, Optional, Callable

logger = logging.getLogger(__name__)

class Topology:
    """
    Defines a signaling network topology for stochastic simulation.
    """

    def __init__(
        self,
        species: List[str],
        parameters: Dict[str, float],
        drift_fn: Optional[Callable] = None,
        name: str = "custom_topology",
        inhibited_species: Optional[str] = None,
        jitted_step_fn: Optional[Callable] = None,
    ):
        self.species = species
        self.parameters = parameters
        self.drift_fn = drift_fn
        self.name = name
        self.inhibited_species = inhibited_species
        self.jitted_step_fn = jitted_step_fn

    @classmethod
    def from_json(cls, json_path: str):
        """Loads a topology from a JSON file."""
        with open(json_path, "r") as f:
            data = json.load(f)
        
        return cls(
            species=data["species"],
            parameters=data["parameters"],
            name=data.get("name", "unnamed_topology")
        )

    @classmethod
    def from_sbml(cls, sbml_path: str, inhibited_species: Optional[str] = None):
        """
        Loads a topology from an SBML file.
        Supports both SBML Core and SBML-qual (qualitative models).
        Requires python-libsbml.
        
        Args:
            sbml_path: Path to the SBML file.
            inhibited_species: ID of the species targeted by the drug.
        """
        try:
            import libsbml
        except ImportError:
            raise ImportError(
                "python-libsbml is required for SBML support. "
                "Install it with `pip install python-libsbml`."
            )

        reader = libsbml.SBMLReader()
        document = reader.readSBML(sbml_path)
        if document.getNumErrors() > 0:
            raise ValueError(f"Error reading SBML file: {document.getErrorLog().toString()}")

        model = document.getModel()
        if model is None:
            raise ValueError("No model found in SBML file.")

        # Check for qualitative model (SBML-qual)
        qual_ext = model.getPlugin("qual")
        if qual_ext:
            return cls._from_sbml_qual(model, qual_ext, inhibited_species)

        species = [s.getId() for s in model.getListOfSpecies()]
        
        # Extract parameters
        parameters = {}
        for p in model.getListOfParameters():
            parameters[p.getId()] = p.getValue()
            
        # If inhibited_species is provided, ensure it exists
        if inhibited_species and inhibited_species not in species:
            logger.warning(f"Inhibited species {inhibited_species} not found in SBML model.")
            
        name = model.getName() or model.getId() or "sbml_topology"
        logger.info(f"Imported SBML model: {name}")

        # Attempt to create a drift function using RoadRunner if available
        drift_fn = None
        try:
            import roadrunner
            def rr_drift_fn(state, params, inhibition=0.0, feedback=None):
                rr = roadrunner.RoadRunner(sbml_path)
                rr.reset()
                # Map state to species
                for i, s_id in enumerate(species):
                    rr[s_id] = state[i]
                
                # Apply inhibition if target is known
                if inhibited_species:
                    rr[inhibited_species] *= (1.0 - inhibition)
                
                # Apply feedback if provided (this is tricky for generic SBML)
                # For now, we scale all rates if feedback is a scalar
                if feedback is not None and isinstance(feedback, (float, int)):
                    rr.model.setGlobalParameterValues(rr.model.getGlobalParameterIds(), 
                                                     [v * feedback for v in rr.model.getGlobalParameterValues()])

                # Get rates of change
                return rr.getRatesOfChange()
            
            drift_fn = rr_drift_fn
            logger.info("Successfully formalized SBML logic using libRoadRunner.")
        except ImportError:
            logger.debug("libRoadRunner not found, SBML drift will use default decay logic.")

        return cls(
            species=species,
            parameters=parameters,
            name=name,
            inhibited_species=inhibited_species,
            drift_fn=drift_fn
        )

    @classmethod
    def _from_sbml_qual(cls, model, qual_ext, inhibited_species: Optional[str] = None):
        """Internal helper to parse SBML-qual models."""
        species = [s.getId() for s in qual_ext.getListOfQualitativeSpecies()]
        parameters = {"k_drift": 0.5, "k_decay": 0.1} # Default rates for logical transitions
        
        name = model.getName() or model.getId() or "sbml_qual_topology"
        logger.info(f"Imported SBML-qual model: {name} with {len(species)} species.")
        
        # We store the transitions for potential use in a drift function
        # For now, we return a Topology that will use generic logic or can be extended
        return cls(
            species=species,
            parameters=parameters,
            name=name,
            inhibited_species=inhibited_species
        )

    def get_initial_state(self) -> np.ndarray:
        """Returns a default initial state (all 0.5 for normalized)."""
        return np.full(len(self.species), 0.5)

def get_default_topology() -> Topology:
    """Returns the default PI3K/AKT/mTOR topology."""
    return Topology(
        species=["PI3K", "AKT", "mTOR"],
        parameters={
            "k_pi3k_base": 0.1,
            "k_pi3k_deg": 0.1,
            "k_akt_act": 0.5,
            "k_akt_deact": 0.1,
            "k_mtor_act": 0.5,
            "k_mtor_deact": 0.1,
        },
        name="PI3K_AKT_mTOR"
    )
