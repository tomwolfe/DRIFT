import unittest
import tempfile
import os
from drift.config import SimulationConfig


class TestSimulationConfig(unittest.TestCase):
    def test_default_config(self):
        """Test that default configuration is created correctly."""
        config = SimulationConfig()
        self.assertEqual(config.drug_kd, 0.5)
        self.assertEqual(config.drug_concentration, 2.0)
        self.assertEqual(config.sim_steps, 100)
        self.assertEqual(config.mc_iterations, 30)
        self.assertEqual(config.dt, 0.1)
        self.assertEqual(config.noise_scale, 0.03)
        self.assertEqual(config.model_name, 'textbook')
        self.assertEqual(config.n_jobs, -1)

    def test_config_validation(self):
        """Test that invalid configurations raise errors."""
        with self.assertRaises(ValueError):
            SimulationConfig(drug_kd=0)  # Kd must be positive
        
        with self.assertRaises(ValueError):
            SimulationConfig(drug_concentration=-1)  # Conc must be non-negative
        
        with self.assertRaises(ValueError):
            SimulationConfig(sim_steps=0)  # Steps must be positive
        
        with self.assertRaises(ValueError):
            SimulationConfig(mc_iterations=0)  # Iterations must be positive
        
        with self.assertRaises(ValueError):
            SimulationConfig(dt=0)  # dt must be positive
        
        with self.assertRaises(ValueError):
            SimulationConfig(noise_scale=-1)  # Noise scale must be non-negative

    def test_from_dict(self):
        """Test creating config from dictionary."""
        config_dict = {
            'drug_kd': 1.0,
            'drug_concentration': 3.0,
            'sim_steps': 50,
            'mc_iterations': 20,
            'dt': 0.05,
            'noise_scale': 0.01,
            'model_name': 'custom_model',
            'n_jobs': 2
        }
        config = SimulationConfig.from_dict(config_dict)
        self.assertEqual(config.drug_kd, 1.0)
        self.assertEqual(config.drug_concentration, 3.0)
        self.assertEqual(config.sim_steps, 50)
        self.assertEqual(config.mc_iterations, 20)
        self.assertEqual(config.dt, 0.05)
        self.assertEqual(config.noise_scale, 0.01)
        self.assertEqual(config.model_name, 'custom_model')
        self.assertEqual(config.n_jobs, 2)

    def test_to_dict(self):
        """Test converting config to dictionary."""
        config = SimulationConfig(drug_kd=1.0, drug_concentration=3.0)
        config_dict = config.to_dict()
        self.assertEqual(config_dict['drug_kd'], 1.0)
        self.assertEqual(config_dict['drug_concentration'], 3.0)

    def test_from_json_file(self):
        """Test loading config from JSON file."""
        config_dict = {
            'drug_kd': 1.5,
            'drug_concentration': 2.5,
            'sim_steps': 75
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            import json
            json.dump(config_dict, f)
            temp_file = f.name
        
        try:
            config = SimulationConfig.from_json_file(temp_file)
            self.assertEqual(config.drug_kd, 1.5)
            self.assertEqual(config.drug_concentration, 2.5)
            self.assertEqual(config.sim_steps, 75)
        finally:
            os.unlink(temp_file)


if __name__ == '__main__':
    unittest.main()