from drift.workbench import Workbench
import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    This example shows how to perform a parameter sweep to analyze the 
    dose-response relationship between drug concentration and growth rate.
    """
    print("Starting Dose-Response Parameter Sweep...")

    concentrations = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]
    results = []

    for conc in concentrations:
        print(f"Testing concentration: {conc}")
        wb = Workbench(drug_kd=0.5, drug_concentration=conc)
        
        # Run a small ensemble to get a stable mean
        # n_jobs=1 for simplicity in this example
        mc_results = wb.run_monte_carlo(n_sims=10, steps=100, n_jobs=1)
        
        # Calculate mean growth across all simulations and time steps
        all_growths = [h['growth'] for h in mc_results['histories']]
        mean_growth = np.mean(all_growths)
        results.append(mean_growth)

    # Visualization
    plt.figure(figsize=(8, 5))
    plt.plot(concentrations, results, 'o-', linewidth=2, markersize=8)
    plt.xscale('log') if concentrations[1] > 0 else None
    plt.xlabel('Drug Concentration')
    plt.ylabel('Average Growth Rate')
    plt.title('Dose-Response Curve: Drug Concentration vs. Metabolic Output')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    
    print("\nResults Summary (Concentration -> Growth):")
    for c, r in zip(concentrations, results):
        print(f"  {c:.2f} -> {r:.4f}")

    # plt.show() # Uncomment to see the plot

if __name__ == "__main__":
    main()
