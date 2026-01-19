# DRIFT: Multi-Scale Stochastic Research Workbench

[![CI](https://github.com/tomwolfe/DRIFT/actions/workflows/ci.yml/badge.svg)](https://github.com/tomwolfe/DRIFT/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Coverage Status](https://coveralls.io/repos/github/tomwolfe/DRIFT/badge.svg?branch=main)](https://coveralls.io/github/tomwolfe/DRIFT?branch=main)

DRIFT (**D**rug-target **R**esponse **I**ntegrated **F**lux **T**rajectory) is a multi-scale stochastic framework designed to bridge the gap between molecular binding events and systemic metabolic phenotypes.

## 📐 System Architecture

DRIFT integrates three distinct biological scales into a unified simulation loop:

```mermaid
graph TD
    A[Drug Concentration] -->|Hill Equation| B(BindingEngine)
    B -->|Inhibition Factor| C(StochasticIntegrator)
    C -->|Langevin Dynamics| D{Signaling State}
    D -->|PI3K/AKT/mTOR| E(MetabolicBridge)
    E -->|Vmax Constraints| F(DFBASolver)
    F -->|FBA Optimization| G[Metabolic Fluxes & Growth]
    G -->|Feedback/Recording| H[History/Dashboard]
```

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
- **Scientific Validation:** Benchmarked against known biological responses and numerical stability tests.

## 📚 Documentation

Detailed documentation is available in the `docs/` directory:

| Document | Description |
| --- | --- |
| [**System Architecture**](docs/architecture.md) | Deep dive into the SDE solvers and FBA coupling. |
| [**Getting Started Tutorial**](docs/tutorial.md) | A step-by-step guide for new users. |
| [**Validation & Benchmarks**](docs/validation.md) | Quantitative comparison with published models. |
| [**Scientific Validation**](docs/scientific_validation.md) | Detailed verification of biological principles. |
| [**API Reference**](docs/api_reference.md) | Auto-generated class and method documentation. |
| [**Reproducibility Guide**](docs/reproducibility.md) | Docker, Conda, and environment configuration. |
| [**Tool Comparison**](docs/comparison_tools.md) | How DRIFT compares to DeepPurpose, COBRApy, etc. |

## 🧪 Scientific Foundation

DRIFT's modeling approach is grounded in established systems biology literature:

- **Signaling Dynamics:** The PI3K/AKT/mTOR SDE model is inspired by the topology described in *Chen et al. (2009), "Input-output behavior of ErbB signaling pathways."*
- **Metabolic Modeling:** Flux Balance Analysis (FBA) implementation follows the standard protocols in *Orth et al. (2010), "What is flux balance analysis?"*
- **Multi-scale Coupling:** The integration of signaling and metabolism reflects recent advances in dynamic FBA as reviewed in *Reimers et al. (2021).*

## 🚀 Quick Start

### Installation

**Using Pip:**
```bash
pip install -e .
```

**Using Conda:**
```bash
conda env create -f environment.yml
conda activate drift
```

**Using Docker:**
```bash
docker build -t drift .
docker run -it drift
```

## 🧪 Development and Testing

```bash
# Run the full test suite
pytest tests/

# Run tests with coverage
pytest tests/ --cov=drift --cov-report=html

# Run specific test file
pytest tests/test_workbench.py
```

### Current Test Coverage
- **Core Modules:** 85%+ line coverage
- **Workbench Class:** Full parameter validation testing
- **Simulation Engines:** Boundary condition and error handling tests
- **Visualization:** Dashboard generation and data integrity checks

## 📈 Example Results

Running a typical simulation produces:
- Interactive HTML dashboard with multi-scale visualization
- Quantitative summary of growth inhibition and metabolic drift
- Statistical analysis of ensemble variability
- Parameter sensitivity reports

![Sample Dashboard](assets/sample_output.png)

## 🗺️ Project Roadmap

- [x] **v0.1.0:** Initial engine release & Monte Carlo integration.
- [ ] **v0.2.0:** Support for multi-target inhibition (Polypharmacology).
- [ ] **v0.3.0:** Integration with SBML models for custom signaling topologies.
- [ ] **v1.0.0:** Stable release with full experimental dataset benchmarking.

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details on:
- Setting up your development environment
- Code style guidelines
- Submitting pull requests
- Reporting issues

## 📜 Citing DRIFT

If you use DRIFT in your research, please cite it as:
> Wolfe, T. (2025). DRIFT: A Multi-Scale Stochastic Framework for Predicting Drug-Induced Metabolic Drift. GitHub Repository. https://github.com/tomwolfe/DRIFT

## 🛠️ Technical Stack
- **COBRApy:** Flux Balance Analysis.
- **Numba:** JIT-compiled SDE solvers.
- **Plotly:** Multi-scale visualization.
- **Numpy/Scipy:** Numerical backend.

## 📋 Release Notes

**Latest Release (v0.1.0):**
- Initial public release
- Multi-scale integration engine
- Monte Carlo uncertainty quantification
- Interactive visualization dashboard

For a complete history of changes, see our [Changelog](CHANGELOG.md).

## 🚀 Release Process

For information about creating new releases, see our [Release Process](RELEASE_PROCESS.md) documentation.