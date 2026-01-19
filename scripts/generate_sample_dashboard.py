#!/usr/bin/env python3
"""
Script to generate a sample DRIFT dashboard.
Useful for previewing the workbench results without manual configuration.
"""

import os
import sys

# Add the project root to the path so we can import drift
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from drift.workbench import Workbench
from drift.visualization import create_dashboard
from drift.config import SimulationConfig

def main():
    print("[*] Generating sample DRIFT simulation results...")
    
    # Use default config
    config = SimulationConfig(
        drug_kd=0.5,
        drug_concentration=2.0,
        sim_steps=100,
        mc_iterations=20
    )
    
    # Ensure output directory exists
    os.makedirs("outputs", exist_ok=True)
    
    # Initialize Workbench
    workbench = Workbench(
        drug_kd=config.drug_kd,
        drug_concentration=config.drug_concentration,
        model_name=config.model_name,
    )

    print(f"[*] Running {config.mc_iterations} simulations...")
    results = workbench.run_monte_carlo(
        n_sims=config.mc_iterations, 
        steps=config.sim_steps
    )

    print("[*] Creating dashboard...")
    fig = create_dashboard(results)
    
    output_path = "outputs/sample_dashboard.html"
    fig.write_html(output_path)
    
    print(f"[SUCCESS] Dashboard generated at: {os.path.abspath(output_path)}")
    print("[*] You can open this file in your web browser to view the results.")

if __name__ == "__main__":
    main()
