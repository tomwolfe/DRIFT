# Tutorial: Simulating Drug-Induced Metabolic Drift

This tutorial walks you through a complete research workflow using DRIFT, from setting up a drug-target pair to analyzing the metabolic consequences.

## Prerequisites
Ensure you have DRIFT installed in your environment:
```bash
pip install -e .
```

## Step 1: Defining the Experiment
We want to simulate a drug with a $K_d$ of 0.5 nM acting on a cell population at a concentration of 1.2 nM.

```python
from drift.workbench import Workbench

# Initialize the workbench
wb = Workbench(drug_kd=0.5, drug_concentration=1.2)
```

## Step 2: Running a Single Trajectory
A single trajectory shows how one "cell" responds over time. This is useful for debugging signaling dynamics.

```python
history = wb.run_simulation(steps=200)

# The history object contains:
# - time: array of time steps
# - signaling: [PI3K, AKT, mTOR] trajectories
# - growth: the resulting biomass flux (growth rate)
```

## Step 3: Assessing Uncertainty with Monte Carlo
Because biological systems are stochastic, we need to run an ensemble of simulations to see the "drift."

```python
results = wb.run_monte_carlo(n_sims=50, steps=200)

# 'results' contains a list of all histories and the basal growth rate
print(f"Mean Growth: {sum(h['growth'][-1] for h in results['histories']) / 50}")
```

## Step 4: Visualizing the Results
DRIFT makes it easy to generate a professional dashboard.

```python
from drift.visualization import create_dashboard

# This will generate 'drift_dashboard.html'
create_dashboard(results)
```

## Analyzing the Output
When you open `drift_dashboard.html`, look for:
1. **Signaling Convergence:** Do the protein levels reach a new (lower) steady state?
2. **Growth Envelopes:** The shaded area in the growth plot represents the 95% confidence interval of the metabolic response. A wide envelope indicates high sensitivity to stochastic signaling noise.
3. **Flux Trajectories:** How quickly does the metabolism "drift" from its optimal state to the drug-inhibited state?

## Next Steps
- Try modifying the kinetic parameters in `drift/signaling.py`.
- Swap the metabolic model in the `Workbench` constructor (e.g., `model_name='e_coli_core'`).
