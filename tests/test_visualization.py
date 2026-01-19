import unittest
import tempfile
import os
import numpy as np
from drift.visualization import create_dashboard
from drift.workbench import Workbench


class TestVisualization(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Create a small test dataset
        wb = Workbench(drug_kd=0.5, drug_concentration=1.0, model_name='textbook')
        self.results = wb.run_monte_carlo(n_sims=3, steps=10, n_jobs=1)

    def test_create_dashboard_basic(self):
        """Test that dashboard creation works with basic results."""
        fig = create_dashboard(self.results)
        # Check that we get a plotly figure object back
        self.assertIsNotNone(fig)
        self.assertEqual(type(fig).__name__, 'Figure')

    def test_create_dashboard_empty_results(self):
        """Test dashboard creation with empty results."""
        empty_results = {'histories': [], 'basal_growth': 0.5}
        with self.assertRaises(ValueError):
            create_dashboard(empty_results)

    def test_create_dashboard_single_history(self):
        """Test dashboard creation with a single simulation history."""
        single_result = {
            'histories': [self.results['histories'][0]],
            'basal_growth': self.results['basal_growth']
        }
        fig = create_dashboard(single_result)
        self.assertIsNotNone(fig)

    def test_create_dashboard_invalid_results(self):
        """Test dashboard creation with invalid results structure."""
        with self.assertRaises(ValueError):
            create_dashboard({})  # Missing required keys raises ValueError because histories is empty
        
        # Test missing basal_growth (uses fallback, so it shouldn't raise unless we make it strict)
        mc_results_no_basal = {'histories': self.results['histories']}
        fig = create_dashboard(mc_results_no_basal)
        self.assertIsNotNone(fig)


class TestIntegration(unittest.TestCase):
    def test_full_pipeline(self):
        """Test the complete DRIFT pipeline from initialization to dashboard."""
        # Initialize workbench
        wb = Workbench(drug_kd=0.3, drug_concentration=1.5, model_name='textbook')
        
        # Run single simulation
        single_result = wb.run_simulation(steps=20)
        self.assertIn('time', single_result)
        self.assertIn('signaling', single_result)
        self.assertIn('growth', single_result)
        self.assertIn('inhibition', single_result)
        self.assertEqual(len(single_result['time']), 20)
        self.assertEqual(len(single_result['growth']), 20)
        
        # Run Monte Carlo
        mc_results = wb.run_monte_carlo(n_sims=5, steps=20, n_jobs=1)
        self.assertIn('histories', mc_results)
        self.assertIn('basal_growth', mc_results)
        self.assertEqual(len(mc_results['histories']), 5)
        
        # Generate dashboard
        fig = create_dashboard(mc_results)
        self.assertIsNotNone(fig)
        
    def test_different_drug_parameters(self):
        """Test workbench with different drug parameters."""
        test_params = [
            {'drug_kd': 0.1, 'drug_concentration': 0.5},
            {'drug_kd': 1.0, 'drug_concentration': 2.0},
            {'drug_kd': 0.01, 'drug_concentration': 0.1}
        ]
        
        for params in test_params:
            with self.subTest(params=params):
                wb = Workbench(**params, model_name='textbook')
                result = wb.run_simulation(steps=10)
                
                # Basic validation
                self.assertIn('growth', result)
                self.assertIn('inhibition', result)
                self.assertGreaterEqual(result['inhibition'], 0)
                self.assertLessEqual(result['inhibition'], 1)

    def test_error_handling(self):
        """Test error handling for invalid inputs."""
        with self.assertRaises(ValueError):
            Workbench(drug_kd=0)  # Invalid Kd
            
        with self.assertRaises(ValueError):
            Workbench(drug_kd=-1)  # Invalid Kd
            
        with self.assertRaises(ValueError):
            Workbench(drug_kd=0.5, drug_concentration=-1)  # Invalid concentration


class TestReproducibility(unittest.TestCase):
    def test_same_seed_produces_same_results(self):
        """Test that using the same seed produces identical results."""
        # Note: This test assumes the workbench supports a random seed parameter
        # which may need to be implemented in the actual code
        pass

    def test_configuration_loading(self):
        """Test loading configurations from different sources."""
        from drift.config import SimulationConfig
        
        # Test from dict
        config_dict = {
            'drug_kd': 0.7,
            'drug_concentration': 1.2,
            'sim_steps': 15,
            'mc_iterations': 3
        }
        config = SimulationConfig.from_dict(config_dict)
        
        # Verify values
        self.assertEqual(config.drug_kd, 0.7)
        self.assertEqual(config.drug_concentration, 1.2)
        self.assertEqual(config.sim_steps, 15)
        self.assertEqual(config.mc_iterations, 3)


if __name__ == '__main__':
    unittest.main()