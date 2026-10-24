import numpy as np
import logging
from typing import Dict, Any, Optional, List, TypedDict, Union
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver
from .core.reproducibility import get_rng
from .core.audit import SimulationAuditLog, AuditSeverity
from .core.exceptions import SolverError

logger = logging.getLogger(__name__)


def _primary_inhibition(inhibition: Union[float, Dict[str, float]]) -> float:
    """The representative inhibition of a (possibly multi-target) result."""
    if isinstance(inhibition, dict):
        return max(inhibition.values()) if inhibition else 0.0
    return float(inhibition)


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
    concentration: np.ndarray
    inhibition_series: np.ndarray
    drug_kd: Union[float, Dict[str, float]]
    params: Dict[str, float]
    headless: bool
    species_names: List[str]
    audit_log: List[Dict[str, Any]]

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

    def _prepare_profile(
        self,
        drug_profile: Optional[np.ndarray],
        drug_profile_times: Optional[np.ndarray],
    ) -> tuple:
        """Validate and normalise the concentration profile.

        Returns ``(profile, times)``; both are ``None`` when the run is
        driven by a scalar concentration instead.
        """
        if drug_profile is None:
            if drug_profile_times is not None:
                raise ValueError(
                    "drug_profile_times was given without a drug_profile")
            return None, None

        profile = np.asarray(drug_profile, dtype=np.float64).reshape(-1)
        if profile.size == 0:
            raise ValueError("drug_profile must not be empty")
        if drug_profile_times is None:
            times = np.arange(profile.size, dtype=np.float64) * self.integrator.dt
        else:
            times = np.asarray(drug_profile_times, dtype=np.float64).reshape(-1)
        if times.size != profile.size:
            raise ValueError(
                "drug_profile and drug_profile_times must have the same "
                f"length ({profile.size} != {times.size})")
        if np.any(np.diff(times) < 0):
            raise ValueError("drug_profile_times must be non-decreasing")
        if np.any(profile < 0):
            raise ValueError("drug_profile must be non-negative")
        return profile, times

    def run(
        self,
        steps: int,
        drug_concentration: Optional[float] = None,
        sub_steps: int = 5,
        refine_feedback: bool = False,
        custom_params: Optional[Dict[str, float]] = None,
        seed: Optional[int] = None,
        drug_profile: Optional[np.ndarray] = None,
        drug_profile_times: Optional[np.ndarray] = None
    ) -> SimulationResult:
        """
        Runs the simulation loop.

        The drug can be given in two ways:

        * ``drug_concentration``: a scalar concentration that drives every
          step with the same value;
        * ``drug_profile`` (with ``drug_profile_times``): the
          concentration-TIME PROFILE of the trial -- VeriTrial's PBPK
          output -- sampled at the step's own time.  This is the seam
          between the two engines: the plasma concentration profile is
          what the PD layer sees.

        Exactly one of the two must be given.
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

        profile, profile_times = self._prepare_profile(
            drug_profile, drug_profile_times)
        if profile is None and drug_concentration is None:
            raise ValueError(
                "either drug_concentration or drug_profile must be given")

        # The exposure of every step: the concentration the PD layer sees
        # and the inhibition it derives from it.
        scalar_concentration = (float(drug_concentration)
                                if drug_concentration is not None else 0.0)

        def step_exposure(step: int) -> tuple:
            if profile is None:
                concentration = scalar_concentration
            else:
                concentration = float(np.interp(
                    step * self.integrator.dt, profile_times, profile))
            return concentration, self.binding.calculate_inhibition(concentration)

        # Calculate a representative inhibition for logging/history
        primary_inhibition = _primary_inhibition(
            self.binding.calculate_inhibition(scalar_concentration)
        )

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
            "concentration": np.zeros(steps),
            "inhibition_series": np.zeros(steps),
            "drug_kd": hist_kd,
            "params": custom_params or {},
            "headless": self.solver.headless,
            "species_names": self.integrator.topology.species,
            "audit_log": []  # Initialize empty audit log
        }

        signaling_hist = []
        growth_hist = []
        concentration_hist = []
        inhibition_hist = []

        for step in range(steps):
            # The PK seam: the concentration of this step, and the
            # inhibition DRIFT's Hill binding derives from it.
            concentration, inhibition = step_exposure(step)
            concentration_hist.append(concentration)
            inhibition_hist.append(_primary_inhibition(inhibition))

            if history["cell_death"]:
                # Once dead, stay dead
                signaling_hist.append(state.copy())
                growth_hist.append(0.0)
                history["status"].append("dead")
                continue

            state_old = state.copy()
            feedback_old = feedback.copy()

            # 1. Predictor signaling step
            for _ in range(sub_steps):
                orig_dt = self.integrator.dt
                self.integrator.dt = dt_small
                state = self.integrator.step(state, inhibition, feedback=feedback)
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
                    state = self.integrator.step(state, inhibition, feedback=avg_feedback)
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
        history["concentration"] = np.asarray(
            concentration_hist, dtype=np.float64)
        history["inhibition_series"] = np.asarray(
            inhibition_hist, dtype=np.float64)
        if profile is not None:
            # The peak exposure of the trial: the largest inhibition the
            # profile ever asks for.
            history["inhibition"] = (
                float(history["inhibition_series"].max())
                if history["inhibition_series"].size else 0.0
            )

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

