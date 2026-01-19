# Validation & Benchmarking

DRIFT is designed to be a research-grade framework. Its design is grounded in established biological literature and validated through quantitative benchmarks.

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

## 4. Comprehensive Scientific Validation
For detailed validation tests including theoretical models, biological principles, and numerical methods, see our [Scientific Validation](scientific_validation.md) document.

## 5. Comparison with Existing Tools
For a detailed comparison of DRIFT with other computational tools in systems biology, see our [Comparison with Existing Tools](comparison_tools.md) document.

## 7. Selected References

DRIFT's architecture is built on foundations established in these key publications:

1.  **Signaling Topologies:** Chen, W. W., et al. (2009). "Input-output behavior of ErbB signaling pathways as revealed by a computer model." *Molecular Systems Biology*, 5(1).
2.  **FBA Methodology:** Orth, J. D., Thiele, I., & Palsson, B. Ø. (2010). "What is flux balance analysis?" *Nature Biotechnology*, 28(3), 245-248.
3.  **Metabolic Control:** Reimers, A. M., et al. (2021). "The role of enzyme kinetics in flux balance analysis." *Current Opinion in Systems Biology*.
4.  **Stochastic Integration:** Higham, D. J. (2001). "An Algorithmic Introduction to Numerical Simulation of Stochastic Differential Equations." *SIAM Review*, 43(3), 525-546.

## 9. Experimental Benchmarking: Alpelisib (BYL719) Case Study

To demonstrate DRIFT's utility in a real-world context, we benchmarked the framework's response against published data for **Alpelisib**, a selective PI3Kα inhibitor used in breast cancer treatment.

### Benchmark Parameters
- **Target:** PI3Kα
- **Experimental IC50:** ~5 nM (in biochemical assays)
- **Cellular Response:** Observed AKT phosphorylation decrease within 30-60 minutes.

### DRIFT Simulation vs. Literature
1. **Dose-Response:** By setting `drug_kd = 0.005` (representing 5 nM), DRIFT reproduces the sigmoidal inhibition of the PI3K/AKT axis consistent with the dose-response curves reported by Fritsch et al. (2014).
2. **Kinetics:** DRIFT's `StochasticIntegrator` captures the rapid dephosphorylation of AKT following PI3K inhibition, matching the 1-hour stabilization window observed in cell lines.
3. **Metabolic Shift:** The resulting reduction in glucose uptake simulated by the `DFBASolver` aligns with PET imaging findings showing decreased FDG-glucose uptake in tumors treated with PI3K inhibitors.

For a runnable version of this benchmark, see `examples/drug_comparison.py`.

## 10. How to Validate Your Own Model

To validate a custom model within DRIFT:
1. **Steady-State Check:** Run a simulation with `drug_concentration = 0` and ensure the growth rate matches the expected $T_d$ (doubling time) of your target cell line.
2. **Inhibition Sweep:** Use `examples/parameter_sweep.py` to generate an $IC_{50}$ curve and compare it with experimental cell viability data.
3. **Trajectory Variance:** Use Monte Carlo simulations to ensure the "stochastic drift" (variance in growth) is within the range of observed single-cell heterogeneity.
