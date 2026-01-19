# Reproducibility Guide

To ensure that DRIFT simulations are consistent across different environments, we provide several tools for reproducibility.

## 1. Docker Environment (Recommended)
Using Docker ensures that all system dependencies (C compilers for Numba, BLAS/LAPACK for COBRApy) are identical.

### Building the Image
```bash
docker build -t drift-workbench .
```

### Running a Simulation
```bash
docker run -v $(pwd)/outputs:/app/outputs drift-workbench
```
This will run a default Monte Carlo ensemble and can be configured to output the `drift_dashboard.html` to your local `outputs/` folder.

## 2. Seeded Randomness
DRIFT uses `numpy.random` for its stochastic signaling noise and Monte Carlo perturbations. For reproducible research:
- All simulations in the `Workbench` use a global seed if set via `numpy.random.seed()`.
- Future versions will support a per-instance `random_state`.

## 3. Configuration Files
Instead of using CLI arguments, we recommend using JSON configuration files for large-scale experiments.
```json
{
  "drug_kd": 0.5,
  "drug_concentration": 1.2,
  "mc_iterations": 100,
  "steps": 200,
  "model_name": "textbook"
}
```
Run with:
```bash
python main.py --config my_experiment.json
```

## 4. Versioning
Ensure you note the version of the metabolic model being used. DRIFT currently supports standard COBRA models (JSON/XML). The default `textbook` model is the *E. coli* core model.
