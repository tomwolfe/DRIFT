import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Any, List, Optional, Union
import logging

logger = logging.getLogger(__name__)

class ExternalValidator:
    """
    Framework for importing and benchmarking simulation results against 
    external experimental datasets (e.g., CSVs of Western Blot or Seahorse data).
    """

    def __init__(self, experimental_data_path: Optional[str] = None):
        """
        Initialize the validator.
        
        Args:
            experimental_data_path: Path to a CSV/Excel file containing experimental data.
        """
        self.exp_data = None
        if experimental_data_path:
            self.load_experimental_data(experimental_data_path)

    def load_experimental_data(self, path: str):
        """Loads experimental data from a file."""
        if path.endswith(".csv"):
            self.exp_data = pd.read_csv(path)
        elif path.endswith(".xlsx") or path.endswith(".xls"):
            self.exp_data = pd.read_excel(path)
        else:
            raise ValueError("Unsupported file format. Use CSV or Excel.")
        logger.info(f"Loaded experimental data from {path}")

    def benchmark(
        self, 
        simulation_history: Dict[str, Any], 
        mapping: Dict[str, str],
        time_col: str = "time"
    ) -> Dict[str, float]:
        """
        Compare simulation history with experimental data.
        
        Args:
            simulation_history: The history dict from a simulation run.
            mapping: Mapping from experimental column names to simulation species/keys.
                     e.g., {"pAKT_level": "AKT", "OCR": "growth"}
            time_col: The column name for time in experimental data.
            
        Returns:
            dict: Benchmark metrics (MSE, correlation) for each mapped entity.
        """
        if self.exp_data is None:
            raise ValueError("No experimental data loaded.")

        results = {}
        sim_time = simulation_history["time"]
        
        from scipy.interpolate import CubicSpline

        for exp_col, sim_key in mapping.items():
            if exp_col not in self.exp_data.columns:
                logger.warning(f"Column {exp_col} not found in experimental data.")
                continue
                
            # Improved interpolation: Cubic Spline for biological curvature
            try:
                # Sort experimental data by time for CubicSpline
                sorted_exp = self.exp_data.sort_values(time_col)
                cs = CubicSpline(sorted_exp[time_col], sorted_exp[exp_col], extrapolate=True)
                exp_val = cs(sim_time)
                # Clip to avoid artifacts below zero if biological data shouldn't be negative
                if (sorted_exp[exp_col] >= 0).all():
                    exp_val = np.maximum(0, exp_val)
            except Exception as e:
                logger.warning(f"CubicSpline interpolation failed for {exp_col}: {e}. Falling back to linear.")
                exp_val = np.interp(
                    sim_time, 
                    self.exp_data[time_col], 
                    self.exp_data[exp_col]
                )
            
            # Extract simulation data
            if sim_key in simulation_history:
                sim_val = simulation_history[sim_key]
            elif "signaling" in simulation_history:
                # Find species index
                from .topology import get_default_topology
                # This is a bit of a hack, ideally history should have species_names
                # For now we assume default topology if not provided
                species = ["PI3K", "AKT", "mTOR"] 
                if sim_key in species:
                    idx = species.index(sim_key)
                    sim_val = simulation_history["signaling"][:, idx]
                else:
                    logger.warning(f"Species {sim_key} not found in signaling history.")
                    continue
            else:
                logger.warning(f"Key {sim_key} not found in simulation history.")
                continue

            # Calculate metrics
            mse = np.mean((sim_val - exp_val)**2)
            corr = np.corrcoef(sim_val, exp_val)[0, 1]
            
            results[sim_key] = {
                "mse": mse,
                "correlation": corr
            }
            
        return results

    def plot_comparison(
        self, 
        simulation_history: Dict[str, Any], 
        mapping: Dict[str, str],
        time_col: str = "time",
        output_path: Optional[str] = None
    ):
        """Visualizes the comparison between simulation and experiment."""
        if self.exp_data is None:
            raise ValueError("No experimental data loaded.")
            
        n_plots = len(mapping)
        fig, axes = plt.subplots(n_plots, 1, figsize=(10, 4 * n_plots), sharex=True)
        if n_plots == 1:
            axes = [axes]
            
        sim_time = simulation_history["time"]
        
        for i, (exp_col, sim_key) in enumerate(mapping.items()):
            ax = axes[i]
            
            # Plot experimental data
            ax.scatter(self.exp_data[time_col], self.exp_data[exp_col], color="red", label="Experimental")
            
            # Plot simulation data
            if sim_key == "growth":
                sim_val = simulation_history["growth"]
            else:
                # Assume signaling
                species = ["PI3K", "AKT", "mTOR"]
                if sim_key in species:
                    idx = species.index(sim_key)
                    sim_val = simulation_history["signaling"][:, idx]
                else:
                    continue
                    
            ax.plot(sim_time, sim_val, label="Simulation", color="blue")
            ax.set_ylabel(sim_key)
            ax.legend()
            
        axes[-1].set_xlabel("Time")
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path)
            logger.info(f"Comparison plot saved to {output_path}")
        else:
            plt.show()
