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

def drift_model(name: str):
    """
    Decorator to mark a function as a drift model for a Topology.
    The function should have the signature: fn(state, params_array, feedback=None)
    where params_array includes the kinetic parameters and the drug inhibition at the end.
    """
    def decorator(fn):
        # We also attempt to JIT it immediately if it's not already
        if type(fn).__name__ != "CPUDispatcher":
            try:
                from numba import njit
                fn = njit(fn)
            except ImportError:
                logger.warning("Numba not installed, drift_model will not be jitted.")
        
        fn._drift_model_name = name
        return fn
    return decorator

@drift_model("PI3K_AKT_mTOR")
def pi3k_akt_mtor_drift(state, params, feedback=None):
    """
    Specific drift function for PI3K/AKT/mTOR with metabolic feedback.
    state: [PI3K, AKT, mTOR]
    params: [k_pi3k_base, k_pi3k_deg, k_akt_act, k_akt_deact, k_mtor_act, k_mtor_deact, inhibition]
    """
    if len(state) < 3:
        # Fallback or error for Numba
        return np.zeros_like(state)

    pi3k, akt, mtor = state[0], state[1], state[2]
    (
        k_pi3k_base,
        k_pi3k_deg,
        k_akt_act,
        k_akt_deact,
        k_mtor_act,
        k_mtor_deact,
        inhibition,
    ) = params[:7]

    # Default feedback to 1.0 if not provided
    fb = np.ones(len(state))
    if feedback is not None:
        fb = feedback

    # Drug inhibits PI3K activity/availability
    effective_pi3k = pi3k * (1.0 - inhibition)

    # PI3K dynamics (basal synthesis and degradation)
    # feedback[0] modulates PI3K synthesis
    dpi3k = k_pi3k_base * fb[0] - k_pi3k_deg * pi3k

    # AKT dynamics (activated by PI3K)
    # feedback[1] modulates AKT activation rate
    dakt = k_akt_act * fb[1] * effective_pi3k * (1.0 - akt) - k_akt_deact * akt

    # mTOR dynamics (activated by AKT, modulated by metabolic feedback)
    # Primary metabolic feedback acts on mTOR (index 2)
    effective_mtor_act = k_mtor_act * fb[2]
    dmtor = effective_mtor_act * akt * (1.0 - mtor) - k_mtor_deact * mtor

    res = np.zeros_like(state)
    res[0] = dpi3k
    res[1] = dakt
    res[2] = dmtor
    
    # Any additional species (like AMPK) just have basal decay in this specific model
    if len(state) > 3:
        for i in range(3, len(state)):
            # Generic decay for unspecified species
            res[i] = 0.1 * (0.5 - state[i])

    return res

@drift_model("Complex_Signaling")
def complex_signaling_drift(state, params, feedback=None):
    """
    Enhanced signaling model with PI3K, AKT, mTOR, and AMPK.
    Includes crosstalk and mechanistic feedback loops.
    state: [PI3K, AKT, mTOR, AMPK]
    params: [k_pi3k_base, k_pi3k_deg, k_akt_act, k_akt_deact, k_mtor_act, k_mtor_deact, k_ampk_act, k_ampk_deact, inhibition]
    """
    if len(state) < 4:
        return np.zeros_like(state)

    pi3k, akt, mtor, ampk = state[0], state[1], state[2], state[3]
    (
        k_pi3k_base, k_pi3k_deg,
        k_akt_act, k_akt_deact,
        k_mtor_act, k_mtor_deact,
        k_ampk_act, k_ampk_deact,
        inhibition
    ) = params[:9]

    fb = np.ones(len(state))
    if feedback is not None:
        fb = feedback

    # 1. PI3K dynamics (inhibited by drug AND mTOR feedback loop)
    # mTORC1/S6K feedback: high mTOR inhibits PI3K activation
    mtor_feedback = 1.0 / (1.0 + 5.0 * mtor) 
    dpi3k = k_pi3k_base * fb[0] * mtor_feedback - k_pi3k_deg * pi3k
    effective_pi3k = pi3k * (1.0 - inhibition)

    # 2. AKT dynamics (activated by PI3K)
    dakt = k_akt_act * fb[1] * effective_pi3k * (1.0 - akt) - k_akt_deact * akt

    # 3. mTOR dynamics (activated by AKT, INHIBITED by AMPK)
    # AMPK is a key metabolic sensor that inhibits mTOR
    ampk_inhibition = 1.0 / (1.0 + 10.0 * ampk)
    dmtor = k_mtor_act * fb[2] * akt * ampk_inhibition * (1.0 - mtor) - k_mtor_deact * mtor

    # 4. AMPK dynamics (inhibited by ATP/metabolic state)
    # feedback[3] represents the energy/ATP probe value
    # High energy -> Low AMPK
    energy_state = fb[3]
    dampk = k_ampk_act * (1.0 - energy_state) * (1.0 - ampk) - k_ampk_deact * ampk

    res = np.zeros_like(state)
    res[0], res[1], res[2], res[3] = dpi3k, dakt, dmtor, dampk
    return res

def get_complex_topology() -> Topology:
    """Returns an advanced 4-node topology with feedback and AMPK sensing."""
    return Topology(
        species=["PI3K", "AKT", "mTOR", "AMPK"],
        parameters={
            "k_pi3k_base": 0.1, "k_pi3k_deg": 0.1,
            "k_akt_act": 0.5, "k_akt_deact": 0.1,
            "k_mtor_act": 0.5, "k_mtor_deact": 0.1,
            "k_ampk_act": 0.3, "k_ampk_deact": 0.1
        },
        drift_fn=complex_signaling_drift,
        name="Complex_Sensing"
    )

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
        drift_fn=pi3k_akt_mtor_drift,
        name="PI3K_AKT_mTOR"
    )
