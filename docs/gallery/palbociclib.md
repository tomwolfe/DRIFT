# Palbociclib (Ibrance) Case Study

## Overview
Palbociclib is a selective inhibitor of cyclin-dependent kinases CDK4 and CDK6, blocking cell cycle progression from G1 to S phase.

## DRIFT Modeling Parameters
- **Target:** CDK4/CDK6
- **Mechanism:** Reversible inhibition of the kinase complex.
- **DRIFT Config:**
  - `drug_kd`: ~10 nM
  - `metabolic_bridge`: Direct coupling between signaling state and the Biomass objective function.

## Predicted Response
1. **Signaling:** Minimal impact on upstream PI3K/AKT signaling, but strong suppression of "Cell Cycle" nodes.
2. **Metabolism:** Decoupling of signaling-driven growth from metabolic flux availability.
3. **Phenotype:** Robust G1 arrest and growth rate reduction.

## Industry Impact
This case study demonstrates DRIFT's flexibility in modeling targets beyond the traditional signaling-metabolism axis, extending into cell cycle control.
