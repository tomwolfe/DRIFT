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
        histories = wb.run_monte_carlo(n_sims=5, steps=5, n_jobs=1) # n_jobs=1 for simplicity in test
        self.assertEqual(len(histories), 5)
        self.assertEqual(len(histories[0]['growth']), 5)

if __name__ == '__main__':
    unittest.main()
