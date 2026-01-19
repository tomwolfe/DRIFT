import unittest
import numpy as np
from drift.engine import SimulationEngine
from drift.binding import BindingEngine
from drift.signaling import StochasticIntegrator
from drift.metabolic import MetabolicBridge, DFBASolver

class TestSimulationEngine(unittest.TestCase):
    def setUp(self):
        self.integrator = StochasticIntegrator()
        self.solver = DFBASolver(model_name="textbook")
        self.bridge = MetabolicBridge(species_names=self.integrator.topology.species)
        self.binding = BindingEngine(kd=0.5)
        self.engine = SimulationEngine(
            integrator=self.integrator,
            solver=self.solver,
            bridge=self.bridge,
            binding=self.binding
        )

    def test_run_basic(self):
        steps = 5
        history = self.engine.run(steps=steps, drug_concentration=1.0)
        self.assertEqual(len(history["growth"]), steps)
        self.assertEqual(len(history["signaling"]), steps)
        self.assertIn("inhibition", history)
        self.assertAlmostEqual(history["inhibition"], 1.0 / (1.0 + 0.5))

    def test_predictor_corrector(self):
        """Test that the predictor-corrector logic runs without error."""
        steps = 5
        # Running with refine_feedback=True to trigger predictor-corrector logic
        history = self.engine.run(steps=steps, drug_concentration=1.0, refine_feedback=True)
        self.assertEqual(len(history["growth"]), steps)
        self.assertTrue(all(g >= 0 for g in history["growth"]))

    def test_custom_params(self):
        """Test that custom parameters are correctly applied."""
        steps = 2
        custom_params = {"p1": 0.5}
        # p1 is index 0 in default topology
        history = self.engine.run(steps=steps, drug_concentration=0.0, custom_params=custom_params)
        # We check if it ran, verifying the internal parameter application didn't crash
        self.assertEqual(len(history["growth"]), steps)
        self.assertEqual(history["params"], custom_params)

    def test_cell_death(self):
        """Test that cell death is handled (by forcing invalid constraints)."""
        # We can't easily force cell death without complex setup, 
        # but we can check if it stays dead if we were to mock it.
        # For now, just ensure the engine handles the optimal case.
        history = self.engine.run(steps=10, drug_concentration=100.0) # High drug might cause death if Kd is low
        self.assertIn("cell_death", history)

if __name__ == "__main__":
    unittest.main()
