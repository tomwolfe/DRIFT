import numpy as np
import logging
from typing import Callable, Dict, Any, Optional, List, TypedDict, Union
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver
from .core.reproducibility import get_rng
from .core.audit import SimulationAuditLog, AuditSeverity
from .core.exceptions import SolverError

logger = logging.getLogger(__name__)

class SimulationResult(TypedDict):
    """Schema for simulation history results."""
    time: np.ndarray
    signaling: np.ndarray
    growth: np.ndarray
    status: List[str]
    cell_death: bool
    death_step: Optional[int]
    death_cause: Optional[str]
    inhibition: float
    drug_kd: Union[float, Dict[str, float]]
    params: Dict[str, float]
    headless: bool
    species_names: List[str]
    audit_log: List[Dict[str, Any]]
    concentration_source: str

class SimulationEngine:
    """
    Dedicated engine for running multi-scale simulations.
    Encapsulates the simulation loop and predictor-corrector logic.
    """

    def __init__(
        self,
        integrator: StochasticIntegrator,
        solver: DFBASolver,
        bridge: MetabolicBridge,
        binding: BindingEngine,
        rng: Optional[np.random.Generator] = None
    ):
        self.integrator = integrator
        self.solver = solver
        self.bridge = bridge
        self.binding = binding
        self.rng = rng or get_rng()

        # Sync RNG with integrator
        self.integrator.rng = self.rng

    def run(
        self,
        steps: int,
        drug_concentration: Union[float, Callable[[float], float], np.ndarray],
        sub_steps: int = 5,
        refine_feedback: bool = False,
        custom_params: Optional[Dict[str, float]] = None,
        seed: Optional[int] = None
    ) -> SimulationResult:
        """
        Runs the simulation loop.

        ``drug_concentration`` may be:
          * a float -- the legacy constant-concentration path;
          * a callable ``C(t)`` -- evaluated at every step time (hours);
          * a 1-D ``numpy`` array -- linearly interpolated over the
            simulation time grid (one sample per step, or any length).

        All three paths produce the same schema; the audit log records which
        one was used so a time-varying run is distinguishable from a float
        placeholder.
        """
        if seed is not None:
            self.rng = get_rng(seed)
            self.integrator.rng = self.rng
            # Also seed Numba for this seed
            from .signaling import seed_numba
            seed_numba(seed)

        # Apply custom signaling parameters if provided
        if custom_params:
            for i, (p_name, p_val) in enumerate(self.integrator.topology.parameters.items()):
                if p_name in custom_params:
                    self.integrator.base_params[i] = custom_params[p_name]

        # Resolve the concentration source once, then evaluate it per step.
        if callable(drug_concentration):
            concentration_source = "callable"
            def _concentration_at(t: float) -> float:
                return float(drug_concentration(t))
        elif isinstance(drug_concentration, np.ndarray):
            concentration_source = "array_interp"
            trace = np.asarray(drug_concentration, dtype=float)
            if trace.ndim != 1 or trace.size == 0:
                raise ValueError("drug_concentration array must be 1-D and non-empty")
            trace_times = np.arange(trace.size) * self.integrator.dt
            def _concentration_at(t: float) -> float:
                return float(np.interp(t, trace_times, trace))
        else:
            concentration_source = "constant"
            constant_concentration = float(drug_concentration)
            def _concentration_at(t: float) -> float:
                return constant_concentration

        inhibition = self.binding.calculate_inhibition(_concentration_at(0.0))

        # Calculate a representative inhibition for logging/history
        if isinstance(inhibition, dict):
            primary_inhibition = max(inhibition.values()) if inhibition else 0.0
        else:
            primary_inhibition = inhibition

        if self.solver.headless:
            logger.warning("SimulationEngine: Running in HEADLESS mode. Metabolism is PROXIED (QUALITATIVE RESULTS ONLY).")

        # Initial state from topology
        state = self.integrator.topology.get_initial_state()
        feedback = np.ones(len(self.integrator.topology.species))  # Initial vector feedback

        # Calculate drug_kd for history: if single target, use float; if multi, use dict.
        if len(self.binding.targets) == 1 and "default" in self.binding.targets:
            hist_kd: Union[float, Dict[str, float]] = self.binding.kd
        else:
            hist_kd = self.binding.targets

        # Initialize audit log
        from .core.audit import SimulationAuditLog
        audit_log = SimulationAuditLog()
        audit_log.log_event(
            step=0,
            component="Drug Concentration",
            message=f"Concentration source: {concentration_source}",
            severity=AuditSeverity.INFO,
            metadata={"source": concentration_source},
        )
        
        dt_small = self.integrator.dt / sub_steps

        if sub_steps > 1:
            audit_log.log_step_size_adjustment(
                step=0,
                old_step_size=self.integrator.dt,
                new_step_size=dt_small,
                reason=f"Sub-stepping enabled (sub_steps={sub_steps}) for integration stability."
            )

        history: SimulationResult = {
            "time": np.arange(steps) * self.integrator.dt,
            "signaling": np.zeros((steps, len(self.integrator.topology.species))),
            "growth": np.zeros(steps),
            "status": [],
            "cell_death": False,
            "death_step": None,
            "death_cause": None,
            "inhibition": primary_inhibition,
            "drug_kd": hist_kd,
            "params": custom_params or {},
            "headless": self.solver.headless,
            "species_names": self.integrator.topology.species,
            "audit_log": [],  # Initialize empty audit log
            "concentration_source": concentration_source,
        }

        signaling_hist = []
        growth_hist = []

        for step in range(steps):
            if history["cell_death"]:
                # Once dead, stay dead
                signaling_hist.append(state.copy())
                growth_hist.append(0.0)
                history["status"].append("dead")
                continue

            step_time = float(history["time"][step])
            step_inhibition = self.binding.calculate_inhibition(
                _concentration_at(step_time))

            state_old = state.copy()
            feedback_old = feedback.copy()

            # 1. Predictor signaling step
            for _ in range(sub_steps):
                orig_dt = self.integrator.dt
                self.integrator.dt = dt_small
                state = self.integrator.step(state, step_inhibition, feedback=feedback)
                self.integrator.dt = orig_dt

            # 2. Metabolic mapping (Signaling -> Metabolism)
            constraints, scalings = self.bridge.get_constraints_with_scalings(state)

            # 3. FBA solver
            try:
                fba_result = self.solver.optimize(constraints, scalings=scalings)
                # Log FBA solver decision
                diagnostic = fba_result.get("diagnostic", "None")
                audit_log.log_solver_decision(step, fba_result["status"], fba_result["objective_value"], diagnostic)
            except SolverError as e:
                audit_log.log_event(
                    step=step,
                    component="FBA Solver",
                    message=f"Critical Solver Error: {str(e)}",
                    severity=AuditSeverity.ERROR,
                    metadata={"error": str(e), "status": "failed"}
                )
                fba_result = {
                    "objective_value": 0.0,
                    "fluxes": {},
                    "status": "failed",
                    "diagnostic": str(e)
                }

            if refine_feedback and fba_result["status"] == "optimal":
                # Predictor-Corrector:
                # Use predicted feedback to re-run signaling for better accuracy
                feedback_pred = self.bridge.get_feedback(fba_result["fluxes"])
                avg_feedback = (feedback_old + feedback_pred) / 2.0

                # Reset and correct signaling
                state = state_old.copy()
                for _ in range(sub_steps):
                    orig_dt = self.integrator.dt
                    self.integrator.dt = dt_small
                    state = self.integrator.step(state, step_inhibition, feedback=avg_feedback)
                    self.integrator.dt = orig_dt

                # Re-run FBA with corrected state
                constraints, scalings = self.bridge.get_constraints_with_scalings(state)
                try:
                    fba_result = self.solver.optimize(constraints, scalings=scalings)
                    # Log the refinement FBA decision
                    audit_log.log_solver_decision(
                        step, 
                        fba_result["status"], 
                        fba_result["objective_value"], 
                        fba_result.get("diagnostic", "Refinement step")
                    )
                except SolverError as e:
                    audit_log.log_event(
                        step=step,
                        component="FBA Solver",
                        message=f"Critical Solver Error during refinement: {str(e)}",
                        severity=AuditSeverity.ERROR,
                        metadata={"error": str(e), "status": "failed"}
                    )
                    fba_result = {
                        "objective_value": 0.0,
                        "fluxes": {},
                        "status": "failed",
                        "diagnostic": str(e)
                    }

                # Log the refinement step
                audit_log.log_event(
                    step=step,
                    component="Predictor-Corrector",
                    message="Applied predictor-corrector refinement",
                    severity=AuditSeverity.INFO
                )

            growth = fba_result["objective_value"]
            fluxes = fba_result["fluxes"]
            status = fba_result["status"]

            if status not in ["optimal", "proxied"]:
                history["cell_death"] = True
                history["death_step"] = step
                diag = fba_result.get("diagnostic", "Unknown constraint")
                logger.warning(f"Cell death detected at step {step} (FBA status: {status}). Cause: {diag}")
                history["death_cause"] = diag
                growth = 0.0
                feedback = np.zeros(len(self.integrator.topology.species))  # Total metabolic collapse
            else:
                feedback = self.bridge.get_feedback(fluxes)

            signaling_hist.append(state.copy())
            growth_hist.append(growth)
            history["status"].append(status)

        history["signaling"] = np.array(signaling_hist)
        history["growth"] = np.array(growth_hist)

        # Convert audit log to serializable format
        history["audit_log"] = [
            {
                "timestamp": event.timestamp.isoformat(),
                "step": event.step,
                "component": event.component,
                "message": event.message,
                "severity": event.severity.value,
                "metadata": event.metadata
            }
            for event in audit_log.events
        ]

        return history

