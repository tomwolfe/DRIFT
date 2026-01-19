# Validation & Benchmarking

DRIFT is designed to be a research-grade framework. While it is currently a computational model, its design is grounded in established biological literature and quantitative benchmarks.

## 1. Biological Consistency
The model's behavior has been validated against known phenotypic responses to PI3K/mTOR inhibition:
- **Dose-Response Sensitivity:** The `BindingEngine` correctly propagates logarithmic drug concentrations into sigmoidal growth inhibition curves, matching standard pharmacological assays.
- **Signaling Cascades:** The `StochasticIntegrator` captures the sequential delay in inhibition from PI3K → AKT → mTOR, as observed in time-resolved phosphoproteomics.
- **Metabolic Shift:** The mapping of mTOR activity to Glucose uptake constraints reflects the Warburg effect's reversal upon treatment with rapalogs or mTOR inhibitors.

## 2. Numerical Stability
- **SDE Convergence:** The Langevin solver (Heun's method equivalent) has been tested for convergence at varying $dt$ values.
- **FBA Feasibility:** The `DFBASolver` includes automated fallbacks to ensure that even under extreme inhibition, the metabolic model remains numerically stable and solvable.

## 3. Comparison with Published Models
The signaling topologies and metabolic constraints in DRIFT are inspired by:
- **Signaling:** Chen et al. (2009) "Input-output behavior of ErbB signaling pathways."
- **Metabolic:** Orth et al. (2010) "What is flux balance analysis?" (Biomass objective function validation).

## 4. How to Validate Your Own Model
To validate a custom model within DRIFT:
1. **Steady-State Check:** Run a simulation with `drug_concentration = 0` and ensure the growth rate matches the expected $T_d$ (doubling time) of your target cell line.
2. **Inhibition Sweep:** Use `examples/parameter_sweep.py` to generate an $IC_{50}$ curve and compare it with experimental cell viability data.
3. **Trajectory Variance:** Use Monte Carlo simulations to ensure the "stochastic drift" (variance in growth) is within the range of observed single-cell heterogeneity.
