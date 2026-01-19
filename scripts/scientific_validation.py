import numpy as np
import matplotlib.pyplot as plt
from drift.workbench import Workbench
from drift.validation import ExternalValidator
import os

def run_scientific_validation():
    print("="*60)
    print("DRIFT SCIENTIFIC VALIDATION: Alpelisib Case Study")
    print("="*60)
    
    # 1. Setup Workbench
    # We use the proxy engine by default if no solver is found, 
    # ensuring the "Headless Gap" is bridged.
    wb = Workbench(drug_kd=1.0, drug_concentration=10.0, model_name="textbook")
    
    exp_data_path = "assets/experimental_alpelisib_data.csv"
    if not os.path.exists(exp_data_path):
        print(f"[!] Error: Experimental data not found at {exp_data_path}")
        return

    # 2. Parameter Calibration
    # Addressing "Parameter Sensitivity" weakness
    print("\n[*] Phase 1: Parameter Calibration")
    best_kd, min_mse = wb.calibrate_to_data(
        experimental_data_path=exp_data_path,
        parameter_name="drug_kd",
        param_range=(0.01, 2.0),
        n_points=15,
        mapping={"growth": "growth", "pAKT": "AKT"}
    )
    
    # 3. Validation Visualization
    # Addressing "Validation vs Verification" weakness
    print("\n[*] Phase 2: Scientific Validation Plotting")
    # Run simulation with best calibrated parameters
    # The experimental data goes up to 48 hours
    steps = int(48.0 / wb.signaling.dt)
    history = wb.run_simulation(steps=steps)
    
    validator = ExternalValidator(exp_data_path)
    output_plot = "outputs/scientific_validation_alpelisib.png"
    os.makedirs("outputs", exist_ok=True)
    
    validator.plot_comparison(
        simulation_history=history,
        mapping={"growth": "growth", "pAKT": "AKT"},
        output_path=output_plot
    )
    print(f"[+] Validation plot saved to {output_plot}")
    
    # 4. Global Sensitivity Analysis
    # Proving the model drivers
    print("\n[*] Phase 3: Global Sensitivity Analysis")
    gsa_results = wb.run_global_sensitivity_analysis(n_sims=40, steps=200)
    
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print(f"  - Calibrated Kd: {best_kd:.4f}")
    print(f"  - Final MSE:     {min_mse:.6f}")
    print(f"  - Headless Mode: {wb.solver.headless}")
    print("="*60)

if __name__ == "__main__":
    run_scientific_validation()
