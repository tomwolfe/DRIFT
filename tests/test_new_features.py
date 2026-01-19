import unittest
import os
import pandas as pd
import warnings
from drift.workbench import Workbench
from drift.signaling import StochasticIntegrator

class TestNewFeatures(unittest.TestCase):
    def test_stability_warning(self):
        """Test that stability warnings are issued for large dt or noise."""
        with self.assertWarns(RuntimeWarning):
            # Large dt
            StochasticIntegrator(dt=0.5, noise_scale=0.02)
        
        with self.assertWarns(RuntimeWarning):
            # High noise impact: 0.4 * sqrt(0.1) = 0.4 * 0.316 = 0.126 > 0.1
            StochasticIntegrator(dt=0.1, noise_scale=0.4)

    def test_monte_carlo_generator(self):
        """Test that run_monte_carlo can return a generator."""
        wb = Workbench(drug_kd=0.5, drug_concentration=1.0, model_name="textbook")
        results_gen = wb.run_monte_carlo(n_sims=3, steps=5, n_jobs=1, return_generator=True)
        
        import types
        self.assertIsInstance(results_gen, types.GeneratorType)
        
        histories = list(results_gen)
        self.assertEqual(len(histories), 3)
        self.assertEqual(len(histories[0]["growth"]), 5)

    def test_incremental_export(self):
        """Test that run_monte_carlo can export incrementally to Parquet."""
        export_path = "tests/test_incremental.parquet"
        if os.path.exists(export_path):
            os.remove(export_path)
            
        wb = Workbench(drug_kd=0.5, drug_concentration=1.0, model_name="textbook")
        # When export_to is provided, it returns a generator
        results_gen = wb.run_monte_carlo(n_sims=3, steps=5, n_jobs=1, export_to=export_path)
        
        # Consume the generator to trigger export
        list(results_gen)
        
        self.assertTrue(os.path.exists(export_path))
        
        # Verify Parquet content
        df = pd.read_parquet(export_path)
        # 3 sims * 5 steps = 15 rows
        self.assertEqual(len(df), 15)
        self.assertIn("sim_id", df.columns)
        self.assertIn("growth", df.columns)
        
        if os.path.exists(export_path):
            os.remove(export_path)

if __name__ == "__main__":
    unittest.main()