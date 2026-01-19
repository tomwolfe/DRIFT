from drift.workbench import Workbench
from drift.visualization import create_dashboard
import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    This example demonstrates Monte Carlo analysis to quantify uncertainty
    and stochastic drift in metabolic responses.
    """
    print("Running Monte Carlo Analysis for Uncertainty Quantification...")

    # Parameters
    DRUG_KD = 0.5
    DRUG_CONCENTRATION = 1.0
    N_SIMULATIONS = 50
    STEPS = 150

    # Initialize Workbench
    wb = Workbench(drug_kd=DRUG_KD, drug_concentration=DRUG_CONCENTRATION)

    # Run Monte Carlo ensemble
    print(f"Running {N_SIMULATIONS} simulations...")
    results = wb.run_monte_carlo(n_sims=N_SIMULATIONS, steps=STEPS, n_jobs=1)

    # Extract data
    histories = results['histories']
    basal_growth = results['basal_growth']

    # Calculate statistics across ensemble
    all_growth_trajectories = np.array([h['growth'] for h in histories])
    mean_growth = np.mean(all_growth_trajectories, axis=0)
    std_growth = np.std(all_growth_trajectories, axis=0)
    ci_lower = mean_growth - 1.96 * std_growth  # 95% confidence interval
    ci_upper = mean_growth + 1.96 * std_growth

    # Plot ensemble results
    plt.figure(figsize=(12, 6))

    # Plot individual trajectories with transparency
    for i, growth_traj in enumerate(all_growth_trajectories):
        plt.plot(growth_traj, alpha=0.3, color='lightblue', linewidth=0.8)

    # Plot mean trajectory with confidence interval
    time_points = range(len(mean_growth))
    plt.plot(mean_growth, color='red', linewidth=2, label='Mean Growth')
    plt.fill_between(time_points, ci_lower, ci_upper, color='red', alpha=0.2, label='95% CI')

    plt.axhline(y=basal_growth, color='black', linestyle='--', alpha=0.7, label='Basal Growth')
    plt.xlabel('Time Steps')
    plt.ylabel('Growth Rate (hr⁻¹)')
    plt.title(f'Monte Carlo Analysis: {N_SIMULATIONS} Stochastic Simulations\nDrug Kd={DRUG_KD}, Conc={DRUG_CONCENTRATION}')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Print summary statistics
    final_growth_means = [h['growth'][-1] for h in histories]
    avg_final_growth = np.mean(final_growth_means)
    std_final_growth = np.std(final_growth_means)
    
    print(f"\nMonte Carlo Summary:")
    print(f"Basal Growth Rate: {basal_growth:.4f}")
    print(f"Mean Final Growth: {avg_final_growth:.4f}")
    print(f"Std of Final Growth: {std_final_growth:.4f}")
    print(f"Growth Inhibition: {(1 - avg_final_growth/basal_growth)*100:.1f}%")
    print(f"Stochastic Variance: {std_final_growth:.4f}")

    # Generate interactive dashboard for detailed analysis
    print("\nGenerating interactive dashboard...")
    dashboard_fig = create_dashboard(results)
    dashboard_fig.write_html("monte_carlo_dashboard.html")
    print("Dashboard saved as monte_carlo_dashboard.html")

    # plt.show() # Uncomment to see the plot

if __name__ == "__main__":
    main()