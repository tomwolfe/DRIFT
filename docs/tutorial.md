# Comprehensive DRIFT Tutorial

This tutorial walks you through a complete research workflow using DRIFT, from setting up a drug-target pair to analyzing the metabolic consequences.

## Prerequisites
Ensure you have DRIFT installed in your environment:
```bash
pip install -e .
```

## Understanding the Framework

DRIFT integrates three key biological scales:
1. **Binding**: Drug-target interactions modeled with dissociation constants ($K_d$)
2. **Signaling**: Stochastic dynamics of signaling cascades (PI3K/AKT/mTOR pathway)
3. **Metabolism**: Flux balance analysis of metabolic networks

The framework connects these scales by:
- Using drug concentration and $K_d$ to determine target occupancy
- Propagating this inhibition through the signaling network
- Mapping signaling states to metabolic constraints
- Solving for growth rates and metabolic fluxes

## Step 1: Basic Simulation Setup

Let's start with a simple simulation of a drug with a $K_d$ of 0.5 µM acting at a concentration of 1.2 µM.

```python
from drift.workbench import Workbench

# Initialize the workbench
wb = Workbench(drug_kd=0.5, drug_concentration=1.2)

# Run a single deterministic simulation
history = wb.run_simulation(steps=200)

# The history object contains:
# - time: array of time steps
# - signaling: [PI3K, AKT, mTOR] trajectories
# - growth: the resulting biomass flux (growth rate)
# - inhibition: average inhibition level
```

## Step 2: Exploring Parameter Space

Different drugs have different potencies. Let's compare how varying $K_d$ affects the response:

```python
import matplotlib.pyplot as plt
import numpy as np

# Compare different drug potencies
kds = [0.1, 0.5, 1.0, 2.0]
colors = ['red', 'orange', 'blue', 'purple']

plt.figure(figsize=(12, 6))

for kd, color in zip(kds, colors):
    wb = Workbench(drug_kd=kd, drug_concentration=1.0)
    history = wb.run_simulation(steps=150)

    plt.plot(history['growth'], label=f'Kd = {kd} µM', color=color)

plt.xlabel('Time Steps')
plt.ylabel('Growth Rate (hr⁻¹)')
plt.title('Effect of Drug Potency on Cellular Growth')
plt.legend()
plt.grid(True, alpha=0.3)
# plt.show()  # Uncomment to display
```

## Step 3: Assessing Uncertainty with Monte Carlo

Biological systems are inherently stochastic. Use Monte Carlo simulations to quantify this uncertainty:

```python
from drift.workbench import Workbench
from drift.visualization import create_dashboard

# Initialize the workbench
wb = Workbench(drug_kd=0.5, drug_concentration=1.2)

# Run an ensemble of 50 simulations
results = wb.run_monte_carlo(n_sims=50, steps=200)

# Extract key statistics
histories = results['histories']
basal_growth = results['basal_growth']

# Calculate ensemble statistics
all_growth_trajectories = [h['growth'] for h in histories]
mean_growth = np.mean(all_growth_trajectories, axis=0)
std_growth = np.std(all_growth_trajectories, axis=0)

print(f"Basal Growth Rate: {basal_growth:.4f}")
print(f"Mean Final Growth: {np.mean([h['growth'][-1] for h in histories]):.4f}")
print(f"Stochastic Variance: {np.std([h['growth'][-1] for h in histories]):.4f}")
```

## Step 4: Advanced Analysis - Dose-Response Relationships

Let's explore how different drug concentrations affect the response:

```python
import numpy as np
import matplotlib.pyplot as plt

def dose_response_analysis():
    concentrations = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    growth_rates = []

    for conc in concentrations:
        wb = Workbench(drug_kd=0.5, drug_concentration=conc)

        # Run Monte Carlo to get stable estimate
        mc_results = wb.run_monte_carlo(n_sims=20, steps=100, n_jobs=1)

        # Calculate mean final growth across all simulations
        final_growths = [h['growth'][-1] for h in mc_results['histories']]
        mean_final_growth = np.mean(final_growths)
        growth_rates.append(mean_final_growth)

    # Plot dose-response curve
    plt.figure(figsize=(10, 6))
    plt.semilogx(concentrations, growth_rates, 'o-', linewidth=2, markersize=8)
    plt.xlabel('Drug Concentration (µM)')
    plt.ylabel('Mean Growth Rate (hr⁻¹)')
    plt.title('Dose-Response Relationship')
    plt.grid(True, which="both", ls="-", alpha=0.5)

    # Calculate approximate IC50
    max_growth = growth_rates[0]  # Growth without drug
    min_growth = min(growth_rates)
    half_effect = max_growth - (max_growth - min_growth) / 2

    # Find concentration closest to half-maximal effect
    closest_idx = min(range(len(growth_rates)),
                     key=lambda i: abs(growth_rates[i] - half_effect))
    approx_ic50 = concentrations[closest_idx]

    print(f"Approximate IC50: {approx_ic50} µM")
    # plt.show()  # Uncomment to display

dose_response_analysis()
```

## Step 5: Visualizing Results with Interactive Dashboards

Generate comprehensive visualizations for detailed analysis:

```python
from drift.visualization import create_dashboard

# Using the results from our Monte Carlo simulation above
dashboard_fig = create_dashboard(results)

# Save the dashboard
dashboard_fig.write_html("comprehensive_analysis.html")
print("Interactive dashboard saved as 'comprehensive_analysis.html'")
```

## Step 6: Custom Configuration

For more control over simulation parameters, use the `SimulationConfig` class:

```python
from drift.config import SimulationConfig

# Create custom configuration
config = SimulationConfig(
    drug_kd=0.3,
    drug_concentration=1.5,
    sim_steps=300,
    mc_iterations=100,
    dt=0.05,  # Smaller time steps for more precision
    noise_scale=0.05,  # Higher noise to emphasize stochastic effects
    model_name='textbook',
    n_jobs=4  # Use 4 parallel processes
)

# Use the configuration to initialize Workbench parameters
wb = Workbench(
    drug_kd=config.drug_kd,
    drug_concentration=config.drug_concentration,
    model_name=config.model_name
)

# Run simulation with custom parameters
results = wb.run_monte_carlo(
    n_sims=config.mc_iterations,
    steps=config.sim_steps,
    dt=config.dt,
    noise_scale=config.noise_scale
)
```

## Analyzing the Output

When you open your dashboard (`comprehensive_analysis.html`), look for:

1. **Signaling Dynamics**: How quickly do the signaling molecules (PI3K, AKT, mTOR) respond to drug treatment? Do they reach a new steady state?

2. **Growth Trajectories**: How does the growth rate change over time? Is there an initial rapid decline followed by stabilization?

3. **Uncertainty Quantification**: The shaded areas represent confidence intervals. Wider intervals indicate greater sensitivity to stochastic effects.

4. **Temporal Patterns**: Are there oscillations or other complex temporal dynamics?

## Troubleshooting Common Issues

### Slow Performance
- Reduce the number of Monte Carlo simulations (`n_sims`)
- Decrease the number of time steps (`steps`)
- Use fewer parallel jobs (`n_jobs`)

### Numerical Instability
- Reduce the time step size (`dt`)
- Lower the noise scale (`noise_scale`)
- Check that your drug parameters are physiologically reasonable

### Memory Issues
- Process simulation results in batches
- Use generators instead of storing all results in memory
- Consider using `n_jobs=1` to limit memory usage

## Advanced Usage Tips

### Custom Metabolic Models
You can use different metabolic models by specifying the `model_name` parameter:

```python
# Use a different model (if available)
wb = Workbench(drug_kd=0.5, drug_concentration=1.0, model_name='e_coli_core')
```

### Parallel Processing
For large ensembles, leverage parallel processing:

```python
# Use all available CPU cores
results = wb.run_monte_carlo(n_sims=100, steps=200, n_jobs=-1)
```

## Next Steps

- Explore the example scripts in the `examples/` directory
- Modify the signaling parameters in `drift/signaling.py` to model different pathways
- Experiment with different metabolic models by changing the `model_name` parameter
- Contribute to the project by adding new features or validation studies
