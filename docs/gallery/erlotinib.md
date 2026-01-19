# Erlotinib (Tarceva) Case Study

## Overview
Erlotinib is a potent epidermal growth factor receptor (EGFR) tyrosine kinase inhibitor used primarily in the treatment of non-small cell lung cancer (NSCLC).

## DRIFT Modeling Parameters
- **Target:** EGFR
- **Mechanism:** Competitive inhibition of the ATP-binding site.
- **DRIFT Config:**
  - `drug_kd`: ~0.5 nM (adjusted for cellular context)
  - `signaling_nodes`: High sensitivity in the PI3K/AKT branch.

## Predicted Response
1. **Signaling:** Rapid suppression of EGFR autophosphorylation and subsequent AKT activation.
2. **Metabolism:** Significant reduction in glycolytic flux, reflecting the dependency of EGFR-mutant cells on aerobic glycolysis.
3. **Phenotype:** Dose-dependent growth arrest (G1 phase), captured by DRIFT's Biomass objective function.

## Visualization
*Simulated trajectories for Erlotinib show high stability in the signaling engine with low stochastic drift compared to multi-target inhibitors.*
