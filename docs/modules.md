# Module Documentation

## `drift.binding`
- `BindingEngine(kd, hill_coefficient)`: Core class for ligand-protein interaction.
  - `calculate_occupancy(conc)`: Returns 0.0-1.0 binding fraction.

## `drift.signaling`
- `StochasticIntegrator(dt, noise_scale)`: High-performance SDE solver.
  - Uses `@njit` (Numba) for the `langevin_step` to ensure Monte Carlo simulations are computationally feasible.
  - Models normalized concentrations of PI3K, AKT, and mTOR.

## `drift.metabolic`
- `MetabolicBridge(reaction_map)`: The interface between signaling and flux.
  - Maps `mTOR` $\to$ `EX_glc__D_e`.
- `DFBASolver(model_name)`: Wrapper for COBRApy.
  - Handles the context-managed optimization of metabolic models.

## `drift.workbench`
- `Workbench`: The orchestrator class.
  - `run_simulation()`: Executes a single time-series.
  - `run_monte_carlo()`: Executes an ensemble for sensitivity analysis.

## `drift.visualization`
- `create_dashboard()`: Generates a 3-panel Plotly figure with:
  - Line plots for signaling.
  - Line plots for flux.
  - Shaded area plots for ensemble statistics.
