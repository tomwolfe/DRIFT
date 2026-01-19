# Cetuximab (Erbitux) Case Study

## Overview
Cetuximab is a monoclonal antibody that binds to the extracellular domain of EGFR, preventing ligand binding and receptor activation.

## DRIFT Modeling Parameters
- **Target:** Extracellular EGFR
- **Mechanism:** Competitive inhibition of ligand (EGF) binding.
- **DRIFT Config:**
  - `drug_kd`: ~0.1-0.5 nM (high affinity)
  - `engine`: BindingEngine configured for antibody kinetics.

## Predicted Response
1. **Signaling:** Inhibition of ligand-induced EGFR dimerization and downstream signaling.
2. **Metabolism:** Reduction in biosynthetic precursors required for cell proliferation.
3. **Phenotype:** Effective growth inhibition in "ligand-driven" virtual cell models.

## Multi-Scale Insight
DRIFT captures the difference between small-molecule TKIs and monoclonal antibodies by modeling the extracellular binding interface and its propagation to the internal metabolic state.
