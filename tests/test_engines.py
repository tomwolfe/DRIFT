import unittest
import numpy as np
import logging
from drift.binding import BindingEngine
from drift.signaling import StochasticIntegrator
from drift.metabolic import MetabolicBridge, DFBASolver

# Suppress noisy INFO logs from external libraries
logging.getLogger("cobra").setLevel(logging.WARNING)
logging.getLogger("httpx").setLevel(logging.WARNING)


class TestBindingEngine(unittest.TestCase):
    def test_occupancy(self):
        be = BindingEngine(kd=0.5)
        # At drug conc = Kd, occupancy should be 0.5
        self.assertAlmostEqual(be.calculate_occupancy(0.5), 0.5)
        self.assertEqual(be.calculate_occupancy(0), 0)
        self.assertGreater(be.calculate_occupancy(1.0), 0.5)

    def test_invalid_kd(self):
        """Test that invalid Kd values raise errors."""
        with self.assertRaises(ValueError):
            BindingEngine(kd=0)
        with self.assertRaises(ValueError):
            BindingEngine(kd=-1)

    def test_invalid_hill_coefficient(self):
        """Test that invalid hill coefficients raise errors."""
        with self.assertRaises(ValueError):
            BindingEngine(kd=1.0, hill_coefficient=0)
        with self.assertRaises(ValueError):
            BindingEngine(kd=1.0, hill_coefficient=-1)

    def test_invalid_concentration(self):
        """Test that invalid concentrations raise errors."""
        be = BindingEngine(kd=0.5)
        with self.assertRaises(ValueError):
            be.calculate_occupancy(-1)

    def test_calculate_inhibition(self):
        """Test that calculate_inhibition works correctly."""
        be = BindingEngine(kd=0.5)
        self.assertEqual(be.calculate_inhibition(0), 0)
        self.assertAlmostEqual(
            be.calculate_inhibition(0.5), 0.5
        )  # Kd = conc => 50% inhibition


class TestSignaling(unittest.TestCase):
    def test_step_shape(self):
        si = StochasticIntegrator()
        initial_state = np.array([0.8, 0.8, 0.8])
        new_state = si.step(initial_state, inhibition=0.5)
        self.assertEqual(new_state.shape, (3,))
        self.assertTrue(np.all(new_state >= 0))
        self.assertTrue(np.all(new_state <= 1))

    def test_invalid_dt_and_noise(self):
        """Test that invalid dt and noise values raise errors."""
        with self.assertRaises(ValueError):
            StochasticIntegrator(dt=0)
        with self.assertRaises(ValueError):
            StochasticIntegrator(dt=-1)
        with self.assertRaises(ValueError):
            StochasticIntegrator(noise_scale=-1)

    def test_invalid_state(self):
        """Test that invalid state raises errors."""
        si = StochasticIntegrator()
        with self.assertRaises(ValueError):
            si.step([0.8, 0.8], 0.5)  # Wrong shape
        with self.assertRaises(ValueError):
            si.step(np.array([0.8]), 0.5)  # Wrong shape
        with self.assertRaises(ValueError):
            si.step("invalid", 0.5)  # Wrong type

    def test_invalid_inhibition(self):
        """Test that invalid inhibition values raise errors."""
        si = StochasticIntegrator()
        with self.assertRaises(ValueError):
            si.step(np.array([0.8, 0.8, 0.8]), -0.1)  # Below 0
        with self.assertRaises(ValueError):
            si.step(np.array([0.8, 0.8, 0.8]), 1.1)  # Above 1


class TestMetabolic(unittest.TestCase):
    def test_bridge_defaults(self):
        bridge = MetabolicBridge()
        # mTOR is index 2. If it's 1.0, scaling should be 1.0 -> -10.0
        constraints = bridge.get_constraints([0.5, 0.5, 1.0])
        self.assertEqual(constraints["EX_glc__D_e"], -10.0)

        # If mTOR is 0.0, scaling should be 0.1 -> -1.0
        constraints = bridge.get_constraints([0.5, 0.5, 0.0])
        self.assertEqual(constraints["EX_glc__D_e"], -1.0)

    def test_bridge_invalid_signaling_state(self):
        """Test that invalid signaling states are handled."""
        bridge = MetabolicBridge()
        # Should NOT raise ValueError for length 2 anymore, just warns
        constraints = bridge.get_constraints([0.5, 0.5])
        self.assertEqual(len(constraints), 0)
        
        with self.assertRaises(ValueError):
            bridge.get_constraints("invalid")  # Wrong type

    def test_bridge_invalid_protein_index(self):
        """Test that invalid protein indices are handled."""
        bridge = MetabolicBridge(
            mappings=[
                {
                    "protein_idx": 5,
                    "reaction_id": "EX_glc__D_e",
                    "influence": "positive",
                    "base_vmax": 10.0,
                }
            ]
        )
        # Should handle invalid index gracefully
        constraints = bridge.get_constraints([0.5, 0.5, 0.5])
        # Since index 5 is invalid, no constraints should be added
        self.assertEqual(len(constraints), 0)

    def test_solver_fallback(self):
        # This might take a moment to download/load models
        with self.assertLogs("drift.metabolic", level="WARNING") as cm:
            solver = DFBASolver(model_name="non_existent_model")
            self.assertEqual(solver.model.id, "e_coli_core")
        self.assertTrue(
            any("Failed to load model 'non_existent_model'" in msg for msg in cm.output)
        )

    def test_solver_optimize(self):
        solver = DFBASolver(model_name="textbook")
        result = solver.solve_step({"EX_glc__D_e": -10.0})
        self.assertGreater(result["objective_value"], 0)
        self.assertIn("EX_glc__D_e", result["fluxes"])

    def test_solver_invalid_constraints(self):
        """Test that invalid constraints raise errors."""
        solver = DFBASolver(model_name="textbook")
        with self.assertRaises(ValueError):
            solver.solve_step("invalid")  # Wrong type


if __name__ == "__main__":
    unittest.main()
