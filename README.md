# DRIFT: Multi-Scale Stochastic Research Workbench

[![CI](https://github.com/tomwolfe/DRIFT/actions/workflows/ci.yml/badge.svg)](https://github.com/tomwolfe/DRIFT/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

DRIFT (**D**rug-target **R**esponse **I**ntegrated **F**lux **T**rajectory) is a multi-scale stochastic framework designed to bridge the gap between molecular binding events and systemic metabolic phenotypes. It enables researchers to simulate how drug-induced signaling perturbations propagate through cellular networks to manifest as metabolic drift.

## 🌟 Key Features

- **Multi-Scale Integration:** Seamlessly couples molecular binding, stochastic signaling (SDEs), and dynamic flux balance analysis (dFBA).
- **Stochastic Dynamics:** Captures cellular heterogeneity using Langevin Dynamics solvers accelerated by Numba.
- **Monte Carlo Uncertainty:** Built-in support for ensemble simulations to assess model robustness and parameter sensitivity.
- **Interactive Dashboards:** Generates comprehensive HTML reports using Plotly for deep-dive analysis of trajectories.
- **Extensible Architecture:** Modular design allows for swapping metabolic models or signaling topologies.

## 🏗️ Architecture

1.  **Molecular Proxy (`BindingEngine`):** Models target occupancy using a probabilistic Hill Equation based on drug concentration and binding affinity ($K_d$).
2.  **Stochastic Signaling (`StochasticIntegrator`):** A high-performance Langevin Dynamics solver that simulates the PI3K/AKT/mTOR signaling axis with intrinsic Gaussian noise.
3.  **Metabolic Phenotype (`DFBASolver`):** A Dynamic Flux Balance Analysis engine that translates signaling states into $V_{max}$ constraints on metabolic reconstructions (e.g., COBRA models).

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/tomwolfe/DRIFT.git
cd DRIFT

# Install as an editable package
pip install -e .
```

### Run Your First Simulation

You can run the workbench via the CLI or use it as a library.

**CLI Usage:**
```bash
python main.py --drug-kd 0.5 --drug-conc 1.0 --mc-iterations 20
```

**Library Usage:**
```python
from drift.workbench import Workbench

# Initialize the workbench
wb = Workbench(drug_kd=0.5, drug_concentration=1.0)

# Run a single trajectory
history = wb.run_simulation(steps=200)

# Access results
print(f"Final Growth Rate: {history['growth'][-1]}")
```

## 📚 Examples & Tutorials

Check out the `examples/` directory for detailed usage:

- [**Basic Simulation**](examples/basic_simulation.py): A walk-through of initializing the workbench and plotting signaling trajectories vs. metabolic output.
- [**Parameter Sweep (Dose-Response)**](examples/parameter_sweep.py): Analyze how varying drug concentrations impact the steady-state growth rate of the model.

## 🔬 Validation

The DRIFT framework has been designed to reflect established biological principles:
- **Signaling Inhibition:** Binding events correctly reduce the activation rate of downstream kinases (AKT, mTOR).
- **Metabolic Coupling:** Reduced mTOR activity is mapped to decreased glucose transporter efficiency and protein synthesis bounds, mimicking the effects of rapalogs.
- **Stochasticity:** The Langevin solver ensures that the metabolic "drift" captures the variance observed in single-cell population studies.

## 🧪 Development and Testing

The project uses `pytest` for unit testing and GitHub Actions for Continuous Integration.

```bash
# Run the test suite
pytest tests/
```

## 🛠️ Technical Stack
- **COBRApy:** Flux Balance Analysis.
- **Numba:** JIT-compiled SDE solvers.
- **Plotly:** Multi-scale visualization.
- **Numpy/Scipy:** Numerical backend.

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.