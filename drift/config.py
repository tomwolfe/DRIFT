"""Configuration module for DRIFT project."""

from dataclasses import dataclass
from typing import Dict, Any
import json


@dataclass
class SimulationConfig:
    """Configuration for DRIFT simulations."""

    # Binding parameters
    drug_kd: float = 0.5  # High affinity binding (lower is stronger)
    drug_concentration: float = 2.0  # Concentration of drug in system

    # Simulation parameters
    sim_steps: int = 100  # Time-steps for simulation
    mc_iterations: int = 30  # Number of Monte Carlo simulations

    # Signaling parameters
    dt: float = 0.1  # Time step for integrator
    noise_scale: float = 0.03  # Noise scale for stochastic integrator

    # Model parameters
    model_name: str = "textbook"  # Metabolic model name

    # Multiprocessing
    n_jobs: int = -1  # Number of parallel jobs (-1 for CPU count)

    def __post_init__(self):
        """Validate configuration parameters."""
        if self.drug_kd <= 0:
            raise ValueError(f"drug_kd must be positive, got {self.drug_kd}")

        if self.drug_concentration < 0:
            raise ValueError(
                f"drug_concentration must be non-negative, got {self.drug_concentration}"
            )

        if self.sim_steps <= 0:
            raise ValueError(f"sim_steps must be positive, got {self.sim_steps}")

        if self.mc_iterations <= 0:
            raise ValueError(
                f"mc_iterations must be positive, got {self.mc_iterations}"
            )

        if self.dt <= 0:
            raise ValueError(f"dt must be positive, got {self.dt}")

        if self.noise_scale < 0:
            raise ValueError(
                f"noise_scale must be non-negative, got {self.noise_scale}"
            )

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]):
        """Create a SimulationConfig from a dictionary."""
        # Filter out keys that aren't in the dataclass
        filtered_dict = {
            k: v for k, v in config_dict.items() if k in cls.__annotations__
        }
        return cls(**filtered_dict)

    @classmethod
    def from_json_file(cls, filepath: str):
        """Load configuration from a JSON file."""
        with open(filepath, "r") as f:
            config_dict = json.load(f)
        return cls.from_dict(config_dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to a dictionary."""
        return {
            "drug_kd": self.drug_kd,
            "drug_concentration": self.drug_concentration,
            "sim_steps": self.sim_steps,
            "mc_iterations": self.mc_iterations,
            "dt": self.dt,
            "noise_scale": self.noise_scale,
            "model_name": self.model_name,
            "n_jobs": self.n_jobs,
        }
