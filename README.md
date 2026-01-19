# DRIFT: Multi-Scale Stochastic Research Workbench

DRIFT (**D**rug-target **R**esponse **I**ntegrated **F**lux **T**rajectory) is a bio-informatics tool designed to bridge the gap between molecular binding events and systemic metabolic phenotypes.

## Architecture

The workbench integrates three distinct scales of biological resolution:

1.  **Molecular Proxy (BindingEngine):** Uses a probabilistic Hill Equation to model target occupancy based on drug concentration and binding affinity ($K_d$).
2.  **Stochastic Signaling (StochasticIntegrator):** A Langevin Dynamics solver (accelerated with Numba) that simulates the PI3K/AKT/mTOR signaling axis. It incorporates Gaussian white noise to represent cellular stochasticity and crowding.
3.  **Metabolic Phenotype (DFBASolver):** A Dynamic Flux Balance Analysis engine that translates signaling states into $V_{max}$ constraints on a human metabolic reconstruction (using the `textbook` core model as a proxy).

## Simulation Workflow

The simulation follows an iterative process:
- **Binding:** Calculate inhibition percentage.
- **Signaling:** Update protein concentrations using SDEs.
- **Bridge:** Map mTOR activity to Glucose Transporter (and other) flux bounds.
- **Metabolism:** Solve for optimal growth (objective function) at each time step.
- **Uncertainty:** A Monte Carlo wrapper runs multiple trajectories to visualize the sensitivity of the system to parameter variations.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

Run the main simulation with default parameters:

```bash
python main.py
```

### Custom Parameters

You can customize the simulation parameters using command-line arguments:

```bash
python main.py --drug-kd 0.75 --drug-conc 1.5 --sim-steps 150 --mc-iterations 50 --model-name textbook
```

### Configuration File

Alternatively, you can use a JSON configuration file:

```bash
python main.py --config config.json
```

Example `config.json`:
```json
{
    "drug_kd": 0.75,
    "drug_concentration": 1.5,
    "sim_steps": 150,
    "mc_iterations": 50,
    "dt": 0.1,
    "noise_scale": 0.025,
    "model_name": "textbook",
    "n_jobs": 4
}
```

This will generate `drift_dashboard.html`, an interactive Plotly dashboard containing:
- **Top Panel:** Stochastic signaling trajectories for PI3K, AKT, and mTOR.
- **Middle Panel:** Temporal growth rate fluctuations.
- **Bottom Panel:** The "Global Homeostatic State" showing the mean phenotypic response with a shaded uncertainty envelope.

## Detailed Documentation

For more in-depth information, please refer to the following:
- [System Architecture](docs/architecture.md): Detailed explanation of the multi-scale coupling logic.
- [Module API](docs/modules.md): Technical breakdown of the Python package structure.

## Technical Stack
- **Python 3.10+**
- **COBRApy:** Constraint-based reconstruction and analysis.
- **Numba:** JIT compilation for high-performance SDE integration.
- **Plotly:** Interactive multi-scale visualization.

## Development and Testing

The project includes a comprehensive unit test suite to ensure the reliability of the multi-scale engines.

Run tests using:
```bash
python -m unittest discover tests
```

### Recent Improvements (Pareto Optimized)
- **Configuration System:** Added flexible configuration system allowing parameter customization via JSON files or command-line arguments.
- **Input Validation:** Comprehensive input validation across all modules to prevent runtime errors.
- **Error Handling:** Improved error handling and logging throughout the application.
- **Documentation:** Enhanced docstrings and documentation for all modules.
- **Test Coverage:** Expanded test suite covering edge cases and error conditions.
- **Performance:** Worker processes in Monte Carlo simulations now reuse model instances and integrators, significantly reducing overhead.
- **Robustness:** Enhanced error handling and logging for metabolic model loading and FBA optimization.
- **Flexibility:** `MetabolicBridge` now supports arbitrary mappings between signaling nodes and metabolic reactions with configurable influence (positive/negative).
