from drift.workbench import Workbench
import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    This example demonstrates comparing different drugs with varying affinities
    and their effects on metabolic drift.
    """
    print("Comparing Multiple Drugs with Different Affinities...")

    # Define different drugs with varying Kd values
    drugs = [
        {'name': 'High Affinity Drug', 'kd': 0.1, 'color': 'red'},
        {'name': 'Medium Affinity Drug', 'kd': 0.5, 'color': 'orange'},
        {'name': 'Low Affinity Drug', 'color': 'blue', 'kd': 2.0},
        {'name': 'Very Low Affinity', 'color': 'purple', 'kd': 5.0}
    ]
    
    drug_concentration = 1.0
    steps = 150

    # Store results for each drug
    drug_results = []

    for drug in drugs:
        print(f"Simulating {drug['name']} (Kd = {drug['kd']})...")
        
        # Initialize Workbench with drug parameters
        wb = Workbench(drug_kd=drug['kd'], drug_concentration=drug_concentration)
        
        # Run a single simulation for comparison
        history = wb.run_simulation(steps=steps)
        
        drug_results.append({
            'name': drug['name'],
            'kd': drug['kd'],
            'growth': history['growth'],
            'color': drug['color'],
            'inhibition': history['inhibition']
        })

    # Visualization
    plt.figure(figsize=(12, 8))
    
    # Plot growth trajectories for each drug
    for result in drug_results:
        plt.plot(result['growth'], label=f"{result['name']} (Kd={result['kd']})", 
                color=result['color'], linewidth=2)
    
    plt.xlabel('Time Steps')
    plt.ylabel('Growth Rate (hr⁻¹)')
    plt.title('Comparative Drug Response: Effect of Binding Affinity on Metabolic Drift')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Print comparative summary
    print("\nComparative Summary:")
    print(f"Drug Concentration: {drug_concentration}")
    print("-" * 60)
    print(f"{'Drug Name':<20} {'Kd':<8} {'Final Growth':<12} {'Inhibition %':<12}")
    print("-" * 60)
    
    for result in drug_results:
        final_growth = result['growth'][-1]
        inhibition_pct = result['inhibition'] * 100
        print(f"{result['name'][:19]:<20} {result['kd']:<8.2f} {final_growth:<12.4f} {inhibition_pct:<12.1f}")

    # Additional analysis: IC50-like calculation
    print("\nIC50-like Analysis:")
    growth_rates = [result['growth'][-1] for result in drug_results]
    kds = [result['kd'] for result in drug_results]
    
    # Simple approximation of IC50 from our data points
    # Find the KD where growth is approximately halved from max growth
    max_growth = max(growth_rates)  # Growth without drug (highest KD should be closest)
    min_growth = min(growth_rates)  # Growth with most potent drug
    half_max = (max_growth + min_growth) / 2
    
    # Find closest to half-max growth
    closest_idx = min(range(len(growth_rates)), key=lambda i: abs(growth_rates[i] - half_max))
    approx_ic50 = kds[closest_idx]
    
    print(f"Approximate IC50-like value: {approx_ic50:.2f}")
    
    # plt.show() # Uncomment to see the plot

if __name__ == "__main__":
    main()