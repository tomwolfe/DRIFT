# DRIFT API Reference

This document provides a comprehensive reference for the DRIFT framework's public API.

## Core Classes

### Workbench

The main entry point for running DRIFT simulations.

#### `Workbench(drug_kd, drug_concentration, model_name='textbook', **kwargs)`

Initializes the DRIFT workbench with specified parameters.

**Parameters:**
- `drug_kd` (float): Dissociation constant of the drug (must be positive)
- `drug_concentration` (float): Concentration of the drug (must be non-negative)
- `model_name` (str): Name of the metabolic model to use (default: 'textbook')
- `**kwargs`: Additional parameters passed to the underlying engines

**Raises:**
- `ValueError`: If parameters are invalid

#### `Workbench.run_simulation(steps, dt=None, noise_scale=None)`

Runs a single deterministic or stochastic simulation.

**Parameters:**
- `steps` (int): Number of time steps to simulate (must be positive)
- `dt` (float, optional): Time step size
- `noise_scale` (float, optional): Scale of stochastic noise

**Returns:**
- `dict`: Dictionary containing simulation history with keys:
  - `'time'`: Array of time points
  - `'signaling'`: Array of signaling molecule concentrations [PI3K, AKT, mTOR]
  - `'growth'`: Array of growth rates over time
  - `'inhibition'`: Average inhibition level

**Raises:**
- `ValueError`: If parameters are invalid

#### `Workbench.run_monte_carlo(n_sims, steps, n_jobs=-1, dt=None, noise_scale=None)`

Runs an ensemble of simulations to quantify uncertainty.

**Parameters:**
- `n_sims` (int): Number of simulations in the ensemble (must be positive)
- `steps` (int): Number of time steps per simulation (must be positive)
- `n_jobs` (int): Number of parallel jobs (-1 for all CPU cores)
- `dt` (float, optional): Time step size
- `noise_scale` (float, optional): Scale of stochastic noise

**Returns:**
- `dict`: Dictionary containing:
  - `'histories'`: List of history dictionaries from each simulation
  - `'basal_growth'`: Baseline growth rate without drug

**Raises:**
- `ValueError`: If parameters are invalid

### SimulationConfig

Configuration class for simulation parameters.

#### `SimulationConfig(drug_kd=0.5, drug_concentration=2.0, sim_steps=100, mc_iterations=30, dt=0.1, noise_scale=0.03, model_name='textbook', n_jobs=-1)`

Creates a configuration object with simulation parameters.

**Parameters:**
- `drug_kd` (float): Dissociation constant of the drug
- `drug_concentration` (float): Concentration of the drug
- `sim_steps` (int): Number of simulation steps
- `mc_iterations` (int): Number of Monte Carlo iterations
- `dt` (float): Time step size
- `noise_scale` (float): Scale of stochastic noise
- `model_name` (str): Name of the metabolic model
- `n_jobs` (int): Number of parallel jobs

**Raises:**
- `ValueError`: If parameters are invalid

#### `SimulationConfig.from_dict(config_dict)`

Creates a configuration from a dictionary.

**Parameters:**
- `config_dict` (dict): Dictionary containing configuration parameters

**Returns:**
- `SimulationConfig`: New configuration object

#### `SimulationConfig.to_dict()`

Converts the configuration to a dictionary.

**Returns:**
- `dict`: Dictionary representation of the configuration

#### `SimulationConfig.from_json_file(filepath)`

Loads configuration from a JSON file.

**Parameters:**
- `filepath` (str): Path to the JSON configuration file

**Returns:**
- `SimulationConfig`: New configuration object

## Engine Classes

### BindingEngine

Handles drug-target binding calculations.

#### `BindingEngine(kd, hill_coefficient=1.0)`

Initializes the binding engine.

**Parameters:**
- `kd` (float): Dissociation constant (must be positive)
- `hill_coefficient` (float): Hill coefficient (must be positive, default: 1.0)

**Raises:**
- `ValueError`: If parameters are invalid

#### `BindingEngine.calculate_occupancy(drug_concentration)`

Calculates the fraction of bound targets.

**Parameters:**
- `drug_concentration` (float): Concentration of the drug (must be non-negative)

**Returns:**
- `float`: Fraction of occupied targets [0, 1]

#### `BindingEngine.calculate_inhibition(drug_concentration)`

Calculates the inhibition level based on occupancy.

**Parameters:**
- `drug_concentration` (float): Concentration of the drug (must be non-negative)

**Returns:**
- `float`: Inhibition level [0, 1]

### StochasticIntegrator

Handles stochastic integration of signaling dynamics.

#### `StochasticIntegrator(dt=0.01, noise_scale=0.03)`

Initializes the stochastic integrator.

**Parameters:**
- `dt` (float): Time step size (must be positive, default: 0.01)
- `noise_scale` (float): Scale of stochastic noise (must be non-negative, default: 0.03)

**Raises:**
- `ValueError`: If parameters are invalid

#### `StochasticIntegrator.step(state, inhibition)`

Performs one integration step.

**Parameters:**
- `state` (array-like): Current state [PI3K, AKT, mTOR] (length 3, values in [0,1])
- `inhibition` (float): Inhibition level [0, 1]

**Returns:**
- `numpy.ndarray`: New state after integration step

**Raises:**
- `ValueError`: If parameters are invalid

### SimulationEngine

Runs the multi-scale loop: signaling sub-steps -> metabolic bridge -> FBA,
with the predictor-corrector refinement.

#### `SimulationEngine.run(steps, drug_concentration=None, sub_steps=5, refine_feedback=False, custom_params=None, seed=None, drug_profile=None, drug_profile_times=None)`

Runs the simulation loop.

**Parameters:**
- `steps` (int): Number of simulation steps
- `drug_concentration` (float, optional): Concentration of the drug, applied to every step
- `sub_steps` (int): Signaling sub-steps per simulation step
- `refine_feedback` (bool): Re-run signaling with the metabolic feedback of the solved fluxes
- `custom_params` (dict, optional): Overrides for the signaling parameters
- `seed` (int, optional): Seed for the integrator and the Numba generator
- `drug_profile` (array-like, optional): The plasma concentration PROFILE of the trial (for example VeriTrial's PBPK output).  When given, each step reads the concentration at its own time and derives the inhibition from it -- this is the seam between a PK engine and this PD/metabolic engine.  Without it, `drug_profile_times` defaults to the engine's own time grid.
- `drug_profile_times` (array-like, optional): The times (h) of `drug_profile`; must have the same length and be non-decreasing

Exactly one of `drug_concentration` and `drug_profile` must be given.

**Returns:**
- `dict`: The simulation history -- `time`, `signaling`, `growth`, `status`, `cell_death`, `death_step`, `death_cause`, `inhibition` (the peak inhibition of the run), `concentration` (the exposure of every step), `inhibition_series`, `drug_kd`, `params`, `headless`, `species_names`, `audit_log`

**Raises:**
- `ValueError`: If neither a concentration nor a profile is given, or if the profile is malformed

### MetabolicBridge

Maps signaling states to metabolic constraints.

#### `MetabolicBridge(mappings=None)`

Initializes the metabolic bridge.

**Parameters:**
- `mappings` (list, optional): Custom mappings from signaling proteins to metabolic reactions.  Each mapping carries `reaction_id`, `base_vmax` and the scaling keys; its optional `constraint_kind` selects what the bound means:
  - `capacity` (the default): the mapping CAPS the reaction -- the lower bound is `-(base_vmax * scaling)`, a negative bound that limits an uptake.
  - `requirement`: the mapping states a DEMAND the reaction must carry -- the lower bound is `+(base_vmax * scaling)`.  A demand above what the network can produce makes the LP infeasible, and the engine then reports `cell_death` on the solver's own verdict.

  An unknown `constraint_kind` raises `ValueError`: a typo must not turn a demand into a cap.

#### `MetabolicBridge.get_constraints(signaling_state)`

Gets metabolic constraints based on signaling state.

**Parameters:**
- `signaling_state` (array-like): Signaling state [PI3K, AKT, mTOR] (length 3)

**Returns:**
- `dict`: Dictionary mapping reaction IDs to constraint values

**Raises:**
- `ValueError`: If signaling state is invalid

### DFBASolver

Solves dynamic flux balance analysis problems.

#### `DFBASolver(model_name='textbook')`

Initializes the FBA solver with a metabolic model.

**Parameters:**
- `model_name` (str): Name of the metabolic model to load

#### `DFBASolver.solve_step(constraints)`

Solves the FBA problem for a single time step.

**Parameters:**
- `constraints` (dict): Reaction constraints to apply

**Returns:**
- `tuple`: (growth_rate, flux_dictionary)

## Visualization Functions

### `create_dashboard(results)`

Creates an interactive dashboard from simulation results.

**Parameters:**
- `results` (dict): Results from a Monte Carlo simulation containing 'histories' and 'basal_growth'

**Returns:**
- `plotly.graph_objects.Figure`: Interactive dashboard figure

**Raises:**
- `KeyError`: If results dictionary is missing required keys
- `ValueError`: If results are invalid