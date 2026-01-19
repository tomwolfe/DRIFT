import unittest
from drift.workbench import Workbench

class TestWorkbench(unittest.TestCase):
    def test_run_simulation(self):
        wb = Workbench(drug_kd=0.5, drug_concentration=1.0, model_name='textbook')
        history = wb.run_simulation(steps=10)
        self.assertEqual(len(history['growth']), 10)
        self.assertIn('signaling', history)
        self.assertIn('inhibition', history)

    def test_run_monte_carlo(self):
        wb = Workbench(drug_kd=0.5, drug_concentration=1.0, model_name='textbook')
        results = wb.run_monte_carlo(n_sims=5, steps=5, n_jobs=1) # n_jobs=1 for simplicity in test
        histories = results['histories']
        self.assertEqual(len(histories), 5)
        self.assertEqual(len(histories[0]['growth']), 5)
        self.assertIn('basal_growth', results)

    def test_invalid_parameters(self):
        """Test that invalid parameters raise errors."""
        with self.assertRaises(ValueError):
            Workbench(drug_kd=0)  # Kd must be positive
        with self.assertRaises(ValueError):
            Workbench(drug_kd=-1)  # Kd must be positive
        with self.assertRaises(ValueError):
            Workbench(drug_concentration=-1)  # Conc must be non-negative

    def test_invalid_run_simulation(self):
        """Test that invalid simulation parameters raise errors."""
        wb = Workbench(drug_kd=0.5, drug_concentration=1.0, model_name='textbook')
        with self.assertRaises(ValueError):
            wb.run_simulation(steps=0)  # Steps must be positive
        with self.assertRaises(ValueError):
            wb.run_simulation(steps=-1)  # Steps must be positive

    def test_invalid_run_monte_carlo(self):
        """Test that invalid Monte Carlo parameters raise errors."""
        wb = Workbench(drug_kd=0.5, drug_concentration=1.0, model_name='textbook')
        with self.assertRaises(ValueError):
            wb.run_monte_carlo(n_sims=0, steps=5, n_jobs=1)  # Sims must be positive
        with self.assertRaises(ValueError):
            wb.run_monte_carlo(n_sims=5, steps=0, n_jobs=1)  # Steps must be positive
        with self.assertRaises(ValueError):
            wb.run_monte_carlo(n_sims=-1, steps=5, n_jobs=1)  # Sims must be positive
        with self.assertRaises(ValueError):
            wb.run_monte_carlo(n_sims=5, steps=-1, n_jobs=1)  # Steps must be positive

if __name__ == '__main__':
    unittest.main()
