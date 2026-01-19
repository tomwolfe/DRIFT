import numpy as np
import logging
from typing import Dict, Any, Optional, List
from .binding import BindingEngine
from .signaling import StochasticIntegrator
from .metabolic import MetabolicBridge, DFBASolver

logger = logging.getLogger(__name__)

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
        custom_params: Optional[Dict[str, float]] = None
    ) -> Dict[str, Any]:
        """
        Runs the simulation loop.
        """
        # Apply custom signaling parameters if provided
        if custom_params:
            for i, (p_name, p_val) in enumerate(self.integrator.topology.parameters.items()):
                if p_name in custom_params:
                    self.integrator.base_params[i] = custom_params[p_name]

        inhibition = self.binding.calculate_inhibition(drug_concentration)
        
        # Initial state from topology
        state = self.integrator.topology.get_initial_state()
        feedback = np.ones(len(self.integrator.topology.species))  # Initial vector feedback

        history = {
            "time": np.arange(steps) * self.integrator.dt,
            "signaling": [],
            "growth": [],
            "status": [],
            "cell_death": False,
            "death_step": None,
            "inhibition": inhibition,
            "drug_kd": self.binding.kd,
            "params": custom_params or {},
        }

        dt_small = self.integrator.dt / sub_steps

        for step in range(steps):
            if history["cell_death"]:
                # Once dead, stay dead
                history["signaling"].append(state.copy())
                history["growth"].append(0.0)
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
            constraints = self.bridge.get_constraints(state)
            
            # 3. FBA solver
            fba_result = self.solver.solve_step(constraints)
            
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
                constraints = self.bridge.get_constraints(state)
                fba_result = self.solver.solve_step(constraints)

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

            history["signaling"].append(state.copy())
            history["growth"].append(growth)
            history["status"].append(status)

        history["signaling"] = np.array(history["signaling"])
        history["growth"] = np.array(history["growth"])
        return history
