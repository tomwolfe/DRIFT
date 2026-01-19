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

class TestSignaling(unittest.TestCase):
    def test_step_shape(self):
        si = StochasticIntegrator()
        initial_state = np.array([0.8, 0.8, 0.8])
        new_state = si.step(initial_state, inhibition=0.5)
        self.assertEqual(new_state.shape, (3,))
        self.assertTrue(np.all(new_state >= 0))
        self.assertTrue(np.all(new_state <= 1))

class TestMetabolic(unittest.TestCase):
    def test_bridge_defaults(self):
        bridge = MetabolicBridge()
        # mTOR is index 2. If it's 1.0, scaling should be 1.0 -> -10.0
        constraints = bridge.get_constraints([0.5, 0.5, 1.0])
        self.assertEqual(constraints['EX_glc__D_e'], -10.0)
        
        # If mTOR is 0.0, scaling should be 0.1 -> -1.0
        constraints = bridge.get_constraints([0.5, 0.5, 0.0])
        self.assertEqual(constraints['EX_glc__D_e'], -1.0)

    def test_solver_fallback(self):
        # This might take a moment to download/load models
        with self.assertLogs('drift.metabolic', level='WARNING') as cm:
            solver = DFBASolver(model_name='non_existent_model')
            self.assertEqual(solver.model.id, 'e_coli_core')
        self.assertTrue(any("Failed to load model 'non_existent_model'" in msg for msg in cm.output))
        
    def test_solver_optimize(self):
        solver = DFBASolver(model_name='textbook')
        growth, fluxes = solver.solve_step({'EX_glc__D_e': -10.0})
        self.assertGreater(growth, 0)
        self.assertIn('EX_glc__D_e', fluxes)

if __name__ == '__main__':
    unittest.main()
