import numpy as np
import logging
from typing import Dict, Any, Optional, List, TypedDict, Union
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver

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
        binding: BindingEngine
    ):
        self.integrator = integrator
        self.solver = solver
        self.bridge = bridge
        self.binding = binding

    def run(
        self,
        steps: int,
        drug_concentration: float,
        sub_steps: int = 5,
        refine_feedback: bool = False,
        custom_params: Optional[Dict[str, float]] = None,
        seed: Optional[int] = None
    ) -> SimulationResult:
        """
        Runs the simulation loop.
        """
        if seed is not None:
            self.integrator.set_seed(seed)

        # Apply custom signaling parameters if provided
        if custom_params:
            for i, (p_name, p_val) in enumerate(self.integrator.topology.parameters.items()):
                if p_name in custom_params:
                    self.integrator.base_params[i] = custom_params[p_name]

        inhibition = self.binding.calculate_inhibition(drug_concentration)
        
        # Calculate a representative inhibition for logging/history
        if isinstance(inhibition, dict):
            # Use the mean inhibition as the single value for history if needed, 
            # or just the first one. For history["inhibition"], we'll take the max.
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

        history: SimulationResult = {
            "time": np.arange(steps) * self.integrator.dt,
            "signaling": np.array([]), # Placeholder
            "growth": np.array([]),    # Placeholder
            "status": [],
            "cell_death": False,
            "death_step": None,
            "death_cause": None,
            "inhibition": primary_inhibition,
            "drug_kd": hist_kd,
            "params": custom_params or {},
            "headless": self.solver.headless,
            "species_names": self.integrator.topology.species
        }

        signaling_hist = []
        growth_hist = []

        dt_small = self.integrator.dt / sub_steps

        for step in range(steps):
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
            fba_result = self.solver.solve_step(constraints, scalings=scalings)
            
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
                fba_result = self.solver.solve_step(constraints, scalings=scalings)

            growth = fba_result["objective_value"]
            fluxes = fba_result["fluxes"]
            status = fba_result["status"]

            if status != "optimal":
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
        return history
