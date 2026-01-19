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
    ):
        self.species = species
        self.parameters = parameters
        self.drift_fn = drift_fn
        self.name = name

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
    def from_sbml(cls, sbml_path: str):
        """
        Loads a topology from an SBML file.
        Requires python-libsbml.
        
        Note: This is a simplified importer for the v0.3.0 roadmap.
        It extracts species and basic reaction metadata.
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
        species = [s.getId() for s in model.getListOfSpecies()]
        
        # Extract parameters
        parameters = {}
        for p in model.getListOfParameters():
            parameters[p.getId()] = p.getValue()
            
        logger.info(f"Imported SBML model: {model.getName() or model.getId()}")
        
        return cls(
            species=species,
            parameters=parameters,
            name=model.getName() or model.getId()
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
