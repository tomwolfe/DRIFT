# DRIFT: Multi-Scale Stochastic Research Workbench

[![CI](https://github.com/tomwolfe/DRIFT/actions/workflows/ci.yml/badge.svg)](https://github.com/tomwolfe/DRIFT/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

DRIFT (**D**rug-target **R**esponse **I**ntegrated **F**lux **T**rajectory) is a multi-scale stochastic framework designed to bridge the gap between molecular binding events and systemic metabolic phenotypes.

## ❓ Why DRIFT?

In drug discovery, linking a molecular binding event to a systemic outcome (like growth inhibition) is often treated as a "black box." DRIFT provides a transparent, mechanistic bridge by:
1.  **Capturing Temporal Dynamics:** Moving beyond static $IC_{50}$ values to see how responses evolve.
2.  **Accounting for Stochasticity:** Modeling the "drift" in metabolic states caused by intrinsic cellular noise.
3.  **Integrating Scales:** Coupling pharmacokinetics (Binding), pharmacodynamics (Signaling), and phenotype (Metabolism) in a single unified solver.

## 🌟 Key Features

- **Multi-Scale Integration:** Seamlessly couples molecular binding, stochastic signaling (SDEs), and dynamic flux balance analysis (dFBA).
- **Stochastic Dynamics:** Captures cellular heterogeneity using Langevin Dynamics solvers accelerated by **Numba**.
- **Monte Carlo Uncertainty:** Built-in support for ensemble simulations to assess model robustness and parameter sensitivity.
- **Interactive Dashboards:** Generates comprehensive HTML reports using **Plotly** for deep-dive analysis of trajectories.

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

```python
from drift.workbench import Workbench
from drift.visualization import create_dashboard

# 1. Initialize with drug parameters
wb = Workbench(drug_kd=0.5, drug_concentration=1.0)

# 2. Run an ensemble of 20 simulations
results = wb.run_monte_carlo(n_sims=20, steps=100)

# 3. Generate the interactive dashboard
create_dashboard(results)
```

## 📚 Documentation

- [**System Architecture**](docs/architecture.md): Deep dive into the SDE solvers and FBA coupling.
- [**Getting Started Tutorial**](docs/tutorial.md): A step-by-step guide for new users.
- [**API Reference**](docs/modules.md): Detailed module and class documentation.
- [**Validation & Benchmarks**](docs/validation.md): How we ensure scientific rigor.
- [**Reproducibility Guide**](docs/reproducibility.md): Using Docker and configuration files.

## 🔬 Research & Validation

The DRIFT framework is grounded in established biological principles, mapping signaling axis (PI3K/AKT/mTOR) to metabolic constraints (Glucose uptake/Protein synthesis). For more details on the numerical stability and biological consistency, see our [Validation Document](docs/validation.md).

## 🧪 Development and Testing

```bash
# Run the test suite
pytest tests/
```

## 📜 Citing DRIFT

If you use DRIFT in your research, please cite it as:
> Wolfe, T. (2025). DRIFT: A Multi-Scale Stochastic Framework for Predicting Drug-Induced Metabolic Drift. GitHub Repository. https://github.com/tomwolfe/DRIFT

## 🛠️ Technical Stack
- **COBRApy:** Flux Balance Analysis.
- **Numba:** JIT-compiled SDE solvers.
- **Plotly:** Multi-scale visualization.
- **Numpy/Scipy:** Numerical backend.