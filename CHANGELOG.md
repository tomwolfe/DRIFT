# Changelog

All notable changes to the DRIFT project will be documented in this file.

## [0.2.0] - 2026-01-30

### Added
- Harden polypharmacology validation with strict input validation
- Parallel reproducibility fixes for bit-perfect Monte Carlo simulations
- Explicit 'solver_status' field to exported JSON results

### Changed
- Enhanced BindingEngine to enforce strict validation of targets dictionary
- Fixed seed propagation to Numba-jitted functions in worker processes
- Improved export metadata with solver status information

## [0.1.1] - 2026-01-18

### Added
- Dedicated `outputs/` directory for simulation results.
- `scripts/generate_sample_dashboard.py` for reproducible visualization generation.
- Support for custom drift functions in `Topology` and `StochasticIntegrator`.
- Comprehensive `.gitignore` to maintain repository hygiene.

### Changed
- Moved generated HTML dashboards out of the root directory.
- Refactored `StochasticIntegrator` to be more generic and extensible.
- Enhanced `CONTRIBUTING.md` with more detailed development guidelines.

## [0.1.0] - 2024-01-18

### Added
- Initial release of DRIFT framework
- Multi-scale integration of binding, signaling, and metabolic models
- Stochastic dynamics with Langevin equations
- Monte Carlo ensemble simulations
- Interactive visualization dashboard
- Basic testing framework