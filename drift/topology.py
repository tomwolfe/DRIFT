import numpy as np
import json
import logging
from typing import List, Dict, Any, Optional, Callable, Union

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
    def from_json(cls, json_path: str) -> "Topology":
        """Loads a topology from a JSON file."""
        with open(json_path, "r") as f:
            data = json.load(f)
        
        return cls(
            species=data["species"],
            parameters=data["parameters"],
            name=data.get("name", "unnamed_topology")
        )

    @classmethod
    def from_sbml(cls, sbml_path: str, inhibited_species: Optional[str] = None) -> "Topology":
        """
        Loads a topology from an SBML file.
        Supports both SBML Core and SBML-qual (qualitative models).
        Requires python-libsbml.
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
        
        parameters = {}
        for p in model.getListOfParameters():
            parameters[p.getId()] = p.getValue()
            
        if inhibited_species and inhibited_species not in species:
            logger.warning(f"Inhibited species {inhibited_species} not found in SBML model.")
            
        name = model.getName() or model.getId() or "sbml_topology"
        logger.info(f"Imported SBML model: {name}")

        drift_fn = None
        try:
            import roadrunner
            # Persistent RoadRunner instance for this Topology
            # Note: For multiprocessing, this would need to be re-initialized in each worker
            _rr_instance = roadrunner.RoadRunner(sbml_path)
            
            def rr_drift_fn(state: np.ndarray, params: Any, inhibition: Union[float, np.ndarray] = 0.0, feedback: Optional[np.ndarray] = None) -> Any:
                _rr_instance.reset()
                for i, s_id in enumerate(species):
                    try:
                        _rr_instance[s_id] = state[i]
                    except (KeyError, ValueError, RuntimeError):
                        continue
                
                # Handle inhibition: can be scalar (legacy) or vector
                if isinstance(inhibition, (float, int)):
                    if inhibited_species:
                        _rr_instance[inhibited_species] *= (1.0 - inhibition)
                elif isinstance(inhibition, np.ndarray):
                    for i, s_id in enumerate(species):
                        if inhibition[i] > 0:
                            _rr_instance[s_id] *= (1.0 - inhibition[i])

                if feedback is not None:
                    # Apply feedback to all parameters as a global scaling factor
                    p_ids = _rr_instance.model.getGlobalParameterIds()
                    p_vals = _rr_instance.model.getGlobalParameterValues()
                    # Use mean feedback as a global modulator for metabolic coupling
                    fb_scalar = np.mean(feedback)
                    _rr_instance.model.setGlobalParameterValues(p_ids, [v * fb_scalar for v in p_vals])
                
                return _rr_instance.getRatesOfChange()
                
            drift_fn = rr_drift_fn
            logger.info(f"Successfully formalized SBML logic for '{name}' using libRoadRunner.")
        except ImportError:
            logger.info("libRoadRunner not found. SBML drift will use default decay logic unless jitted_step_fn is provided.")

        return cls(
            species=species,
            parameters=parameters,
            name=name,
            inhibited_species=inhibited_species,
            drift_fn=drift_fn
        )

    @classmethod
    def _from_sbml_qual(cls, model: Any, qual_ext: Any, inhibited_species: Optional[str] = None) -> "Topology":
        """Internal helper to parse SBML-qual models and generate a logic-based drift function."""
        try:
            import libsbml
        except ImportError:
            raise ImportError("libsbml is required for SBML-qual support.")

        name = model.getName() or model.getId() or "sbml_qual_topology"
        species_ids = [s.getId() for s in qual_ext.getListOfQualitativeSpecies()]
        parameters = {"k_drift": 0.5, "k_decay": 0.1}
        
        import re
        def sanitize(sid: str) -> str:
            return re.sub(r'[^a-zA-Z0-9_]', '_', sid)
        
        clean_ids = {s: sanitize(s) for s in species_ids}

        lines = [
            "import numpy as np",
            "from numba import njit",
            "",
            "@njit",
            "def jitted_qual_drift(state, params, inhibition=0.0, feedback=None):",
            "    kd = params[0]",
            "    kn = params[1]",
            "    res = np.zeros(len(state))",
            "    fb = np.ones(len(state))",
            "    if feedback is not None: fb = feedback",
            "",
        ]
        
        for idx, s_id in enumerate(species_ids):
            lines.append(f"    {clean_ids[s_id]} = state[{idx}]")
            # Initialize target levels for all species
            lines.append(f"    target_{clean_ids[s_id]} = 0.0")
        lines.append("")

        for trans in qual_ext.getListOfTransitions():
            outputs = [o.getQualitativeSpecies() for o in trans.getListOfOutputs()]
            default_term = trans.getListOfFunctionTerms().getDefaultTerm()
            default_val = float(default_term.getResultLevel())
            
            # Add function terms (if-elif-else)
            terms = []
            for ft in trans.getListOfFunctionTerms():
                if ft == default_term: continue
                math = ft.getMath()
                formula = libsbml.formulaToL3String(math)
                py_formula = formula
                for s_id in sorted(species_ids, key=len, reverse=True):
                    py_formula = re.sub(rf'\b{s_id}\b', clean_ids[s_id], py_formula)
                
                py_formula = py_formula.replace("&&", " and ").replace("||", " or ").replace("!", " not ")
                if "==" not in py_formula and "=" in py_formula:
                    py_formula = py_formula.replace("=", "==")
                terms.append((py_formula, float(ft.getResultLevel())))
            
            # Start logic block
            # Default value is only used if no terms match
            for out_id in outputs:
                lines.append(f"    target_{clean_ids[out_id]} = {default_val}")
            
            if terms:
                for i, (f, l) in enumerate(terms):
                    prefix = "if" if i == 0 else "elif"
                    lines.append(f"    {prefix} {f}:")
                    for out_id in outputs:
                        lines.append(f"        target_{clean_ids[out_id]} = {l}")
            lines.append("")

        for i, s_id in enumerate(species_ids):
            target_expr = f"target_{clean_ids[s_id]}"
            if s_id == inhibited_species:
                target_expr = f"({target_expr} * (1.0 - inhibition))"
            lines.append(f"    res[{i}] = kd * fb[{i}] * ({target_expr} - state[{i}]) - kn * state[{i}]")
        
        lines.append("    return res")
        
        code = "\n".join(lines)
        local_namespace: Dict[str, Any] = {}
        from numba import njit
        try:
            # JIT compilation of internally generated code from SBML qual model
            compiled_code = compile(code, '<string>', 'exec')
            exec(compiled_code, {"np": np, "njit": njit}, local_namespace)  # nosec B102
            drift_fn = local_namespace["jitted_qual_drift"]
            drift_fn._drift_model_name = f"jitted_{name}"
            logger.info(f"Successfully JIT-compiled SBML-qual drift function for {name}")
        except Exception as e:
            logger.warning(f"Failed to JIT SBML-qual drift: {e}. Falling back to interpreted logic.")
            drift_fn = cls._generate_interpreted_qual_drift(species_ids, qual_ext, inhibited_species)

        return cls(
            species=species_ids,
            parameters=parameters,
            name=name,
            inhibited_species=inhibited_species,
            drift_fn=drift_fn
        )

    @staticmethod
    def _generate_interpreted_qual_drift(species: List[str], qual_ext: Any, inhibited_species: Optional[str]) -> Callable:
        """Fallback interpreted drift function generator."""
        transitions: List[Dict[str, Any]] = []
        import libsbml
        for trans in qual_ext.getListOfTransitions():
            outputs = [o.getQualitativeSpecies() for o in trans.getListOfOutputs()]
            default_term = trans.getListOfFunctionTerms().getDefaultTerm()
            default_val = float(default_term.getResultLevel())
            func_terms = []
            for ft in trans.getListOfFunctionTerms():
                if ft == default_term: continue
                formula = libsbml.formulaToL3String(ft.getMath())
                func_terms.append((formula, float(ft.getResultLevel())))
            
            transitions.append({
                "outputs": outputs,
                "default": default_val,
                "terms": func_terms
            })

        def interpreted_drift(state: np.ndarray, params: Any, inhibition: float = 0.0, feedback: Optional[np.ndarray] = None) -> np.ndarray:
            if isinstance(params, dict):
                kd, kn = params.get("k_drift", 0.5), params.get("k_decay", 0.1)
            else:
                kd, kn = params[0], params[1]

            s_map: Dict[str, int] = {s: i for i, s in enumerate(species)}
            target_levels = np.zeros(len(species))
            import re

            for trans in transitions:  # type: ignore
                active_level = trans["default"]
                for formula, level in trans["terms"]:
                    eval_f = formula
                    for s_id, idx in s_map.items():
                        eval_f = re.sub(rf'\b{s_id}\b', str(state[idx]), eval_f)
                    eval_f = eval_f.replace("&&", " and ").replace("||", " or ").replace("!", " not ")
                    try:
                        # eval is used with restricted namespace for SBML expressions from trusted files
                        result = eval(eval_f, {"__builtins__": {}}, {})  # nosec B307
                        if result:
                            active_level = level
                            break
                    except Exception:  # nosec B110
                        pass
                for out_id in trans["outputs"]:  # type: ignore
                    out_idx = s_map.get(out_id)
                    if out_idx is not None:
                        target_levels[int(out_idx)] = active_level  # type: ignore
            
            res = np.zeros_like(state)
            fb = feedback if feedback is not None else np.ones(len(state))
            for i, s_id in enumerate(species):
                target = target_levels[i]
                if s_id == inhibited_species: target *= (1.0 - inhibition)
                res[i] = kd * fb[i] * (target - state[i]) - kn * state[i]
            return res
            
        return interpreted_drift

    def get_initial_state(self) -> np.ndarray:
        """Returns a default initial state (all 0.5 for normalized)."""
        return np.full(len(self.species), 0.5)

def drift_model(name: str) -> Callable:
    """
    Decorator to mark a function as a drift model for a Topology.
    """
    def decorator(fn: Callable) -> Callable:
        if type(fn).__name__ != "CPUDispatcher":
            try:
                from numba import njit
                fn = njit(fn)
            except ImportError:
                logger.warning("Numba not installed, drift_model will not be jitted.")
        fn._drift_model_name = name # type: ignore
        return fn
    return decorator

@drift_model("PI3K_AKT_mTOR")
def pi3k_akt_mtor_drift(state: np.ndarray, params: np.ndarray, feedback: Optional[np.ndarray] = None) -> np.ndarray:
    if len(state) < 3:
        return np.zeros_like(state)
    pi3k, akt, mtor = state[0], state[1], state[2]
    
    # Generic parameter mapping: first N are base params, last M are inhibitions
    # For this model, we expect 6 base params and 3 inhibition params
    k_pi3k_base, k_pi3k_deg = params[0], params[1]
    k_akt_act, k_akt_deact = params[2], params[3]
    k_mtor_act, k_mtor_deact = params[4], params[5]
    
    # Inhibition can be a single scalar (legacy) or a vector matching species
    if len(params) == 7:
        inh_pi3k = params[6]
        inh_akt = 0.0
        inh_mtor = 0.0
    elif len(params) >= 9:
        inh_pi3k, inh_akt, inh_mtor = params[6], params[7], params[8]
    else:
        inh_pi3k = 0.0
        inh_akt = 0.0
        inh_mtor = 0.0

    fb = feedback if feedback is not None else np.ones(len(state))
    
    effective_pi3k = pi3k * (1.0 - inh_pi3k)
    effective_akt = akt * (1.0 - inh_akt)
    
    dpi3k = k_pi3k_base * fb[0] - k_pi3k_deg * pi3k
    dakt = k_akt_act * fb[1] * effective_pi3k * (1.0 - akt) - k_akt_deact * akt
    dmtor = k_mtor_act * fb[2] * effective_akt * (1.0 - mtor) - k_mtor_deact * mtor
    
    res = np.zeros_like(state)
    res[0], res[1], res[2] = dpi3k, dakt, dmtor
    
    if len(state) > 3:
        for i in range(3, len(state)):
            res[i] = 0.1 * (0.5 - state[i])
    return res

@drift_model("Complex_Signaling")
def complex_signaling_drift(state: np.ndarray, params: np.ndarray, feedback: Optional[np.ndarray] = None) -> np.ndarray:
    if len(state) < 4:
        return np.zeros_like(state)
    pi3k, akt, mtor, ampk = state[0], state[1], state[2], state[3]
    
    # Expect 8 base params and 4 inhibition params
    k_pi3k_base, k_pi3k_deg = params[0], params[1]
    k_akt_act, k_akt_deact = params[2], params[3]
    k_mtor_act, k_mtor_deact = params[4], params[5]
    k_ampk_act, k_ampk_deact = params[6], params[7]

    if len(params) == 9:
        inh_pi3k = params[8]
        inh_akt = inh_mtor = inh_ampk = 0.0
    elif len(params) >= 12:
        inh_pi3k, inh_akt, inh_mtor, inh_ampk = params[8], params[9], params[10], params[11]
    else:
        inh_pi3k = inh_akt = inh_mtor = inh_ampk = 0.0

    fb = feedback if feedback is not None else np.ones(len(state))
    mtor_feedback = 1.0 / (1.0 + 5.0 * mtor) 
    dpi3k = k_pi3k_base * fb[0] * mtor_feedback - k_pi3k_deg * pi3k
    
    effective_pi3k = pi3k * (1.0 - inh_pi3k)
    effective_akt = akt * (1.0 - inh_akt)
    
    dakt = k_akt_act * fb[1] * effective_pi3k * (1.0 - akt) - k_akt_deact * akt
    ampk_inhibition = 1.0 / (1.0 + 10.0 * ampk)
    dmtor = k_mtor_act * fb[2] * effective_akt * ampk_inhibition * (1.0 - mtor) - k_mtor_deact * mtor
    
    energy_state = fb[3]
    dampk = k_ampk_act * (1.0 - energy_state) * (1.0 - ampk) - k_ampk_deact * ampk
    
    res = np.zeros_like(state)
    res[0], res[1], res[2], res[3] = dpi3k, dakt, dmtor, dampk
    return res

def get_complex_topology() -> Topology:
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
    return Topology(
        species=["PI3K", "AKT", "mTOR"],
        parameters={
            "k_pi3k_base": 0.1, "k_pi3k_deg": 0.1,
            "k_akt_act": 0.5, "k_akt_deact": 0.1,
            "k_mtor_act": 0.5, "k_mtor_deact": 0.1,
        },
        drift_fn=pi3k_akt_mtor_drift,
        name="PI3K_AKT_mTOR"
    )