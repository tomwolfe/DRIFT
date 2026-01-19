# Module Documentation

This document provides a technical overview of the internal modules within the DRIFT framework.

## `drift.binding`
The `BindingEngine` handles the molecular scale of the simulation.

- **`BindingEngine(kd: float, hill_coefficient: float = 1.0)`**
  - `calculate_occupancy(drug_concentration: float) -> float`: Returns the fraction of target bound by the drug [0, 1] using the Hill Equation.
  - `calculate_inhibition(drug_concentration: float) -> float`: Alias for `calculate_occupancy`, representing the inhibitory effect on signaling.

## `drift.signaling`
The `StochasticIntegrator` manages the temporal signaling dynamics.

- **`StochasticIntegrator(dt: float = 0.1, noise_scale: float = 0.02)`**
  - **Langevin Dynamics:** Solves a system of SDEs for the PI3K/AKT/mTOR axis.
  - **Performance:** Utilizes Numba-accelerated `langevin_step` for high-speed simulation.
  - `step(state: np.ndarray, inhibition: float) -> np.ndarray`: Advances the signaling state by one time-step `dt`.

## `drift.metabolic`
The metabolic modules bridge signaling states to flux distributions.

- **`MetabolicBridge(mappings: list = None)`**
  - Translates signaling protein levels (e.g., mTOR) to metabolic constraints (e.g., Vmax for Glucose uptake).
  - `get_constraints(signaling_state: np.ndarray) -> dict`: Returns a dictionary of reaction IDs and their calculated lower bounds.

- **`DFBASolver(model_name: str = 'textbook')`**
  - A robust wrapper around COBRApy for solving individual FBA steps.
  - Handles model loading and provides fallbacks for common models like `e_coli_core`.
  - `solve_step(constraints: dict) -> (float, dict)`: Solves the model and returns the objective value and flux distribution.

## `drift.workbench`
The `Workbench` is the high-level API for users.

- **`Workbench(drug_kd=1.0, drug_concentration=2.0, model_name='textbook')`**
  - `run_simulation(steps=100) -> dict`: Executes a single integrated multi-scale trajectory.
  - `run_monte_carlo(n_sims=30, steps=100, n_jobs=-1) -> dict`: Runs an ensemble of simulations in parallel, capturing system sensitivity and uncertainty.

## `drift.visualization`
- **`create_dashboard(results: dict)`**
  - Generates a rich, interactive Plotly dashboard.
  - Visualizes signaling trajectories, growth rate fluctuations, and ensemble statistics with shaded uncertainty envelopes.

## Complete API Reference

For detailed API documentation with all parameters, return types, and examples, see the [API Reference](api_reference.md).