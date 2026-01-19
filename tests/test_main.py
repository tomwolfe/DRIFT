import unittest
import sys
import os
from unittest.mock import patch, MagicMock
import tempfile
import json
from main import main, parse_args
from drift.config import SimulationConfig


class TestMainModule(unittest.TestCase):
    def test_parse_args_basic(self):
        """Test argument parsing with basic options."""
        test_args = ["main.py", "--drug-kd", "0.7", "--drug-conc", "1.5"]

        with patch.object(sys, "argv", test_args):
            args = parse_args()
            self.assertEqual(args.drug_kd, 0.7)
            self.assertEqual(args.drug_conc, 1.5)

    def test_parse_args_with_mc_iterations(self):
        """Test argument parsing with Monte Carlo iterations."""
        test_args = ["main.py", "--mc-iterations", "50", "--sim-steps", "200"]

        with patch.object(sys, "argv", test_args):
            args = parse_args()
            self.assertEqual(args.mc_iterations, 50)
            self.assertEqual(args.sim_steps, 200)

    def test_parse_args_with_config_file(self):
        """Test argument parsing with config file."""
        # Create a temporary config file
        config_data = {
            "drug_kd": 0.8,
            "drug_concentration": 1.2,
            "sim_steps": 100,
            "mc_iterations": 25,
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(config_data, f)
            temp_config_file = f.name

        try:
            test_args = ["main.py", "--config", temp_config_file]

            with patch.object(sys, "argv", test_args):
                args = parse_args()
                self.assertEqual(args.config, temp_config_file)
        finally:
            os.unlink(temp_config_file)

    @patch("main.Workbench")
    @patch("main.create_dashboard")
    @patch("builtins.print")
    @patch("webbrowser.open")
    def test_main_execution(
        self, mock_webbrowser, mock_print, mock_create_dashboard, mock_workbench_class
    ):
        """Test main function execution with mocked dependencies."""
        # Mock the workbench instance and its methods
        mock_workbench_instance = MagicMock()
        mock_workbench_instance.run_monte_carlo.return_value = {
            "histories": [{"growth": [0.5, 0.4, 0.3], "inhibition": 0.5}],
            "basal_growth": 0.6,
        }
        mock_workbench_class.return_value = mock_workbench_instance

        # Mock the dashboard creation
        mock_dashboard_fig = MagicMock()
        mock_create_dashboard.return_value = mock_dashboard_fig

        # Create a test config
        config = SimulationConfig(
            drug_kd=0.5, drug_concentration=1.0, mc_iterations=1, sim_steps=3
        )

        # Run main with the config
        main(config)

        # Verify that workbench was initialized with correct parameters
        mock_workbench_class.assert_called_once_with(
            drug_kd=0.5, 
            drug_concentration=1.0, 
            model_name="textbook",
            time_unit="hours",
            concentration_unit="uM"
        )

        # Verify that run_monte_carlo was called with correct parameters
        mock_workbench_instance.run_monte_carlo.assert_called_once_with(
            n_sims=1, steps=3, n_jobs=-1
        )

        # Verify that dashboard was created
        mock_create_dashboard.assert_called_once()

    @patch("main.Workbench")
    @patch("main.create_dashboard")
    @patch("builtins.print")
    @patch("webbrowser.open")
    def test_main_with_exception(
        self, mock_webbrowser, mock_print, mock_create_dashboard, mock_workbench_class
    ):
        """Test main function handles exceptions properly."""
        # Mock the workbench to raise an exception
        mock_workbench_class.side_effect = Exception("Test error")

        config = SimulationConfig()

        # Run main and expect it to re-raise the exception
        with self.assertRaises(Exception) as context:
            main(config)

        self.assertEqual(str(context.exception), "Test error")

    def test_config_creation_from_args(self):
        """Test creating config from command line arguments."""
        test_args = [
            "main.py",
            "--drug-kd",
            "0.6",
            "--drug-conc",
            "1.4",
            "--sim-steps",
            "50",
            "--mc-iterations",
            "10",
            "--model-name",
            "test_model",
        ]

        with patch.object(sys, "argv", test_args):
            args = parse_args()

            # Create config with parsed arguments
            config_kwargs = {}
            if args.drug_kd is not None:
                config_kwargs["drug_kd"] = args.drug_kd
            if args.drug_conc is not None:
                config_kwargs["drug_concentration"] = args.drug_conc
            if args.sim_steps is not None:
                config_kwargs["sim_steps"] = args.sim_steps
            if args.mc_iterations is not None:
                config_kwargs["mc_iterations"] = args.mc_iterations
            if args.model_name is not None:
                config_kwargs["model_name"] = args.model_name

            config = SimulationConfig(**config_kwargs)

            # Verify config values
            self.assertEqual(config.drug_kd, 0.6)
            self.assertEqual(config.drug_concentration, 1.4)
            self.assertEqual(config.sim_steps, 50)
            self.assertEqual(config.mc_iterations, 10)
            self.assertEqual(config.model_name, "test_model")


if __name__ == "__main__":
    unittest.main()
