# Imatinib (Gleevec) Case Study

## Overview
Imatinib revolutionized the treatment of Chronic Myeloid Leukemia (CML) by specifically targeting the BCR-ABL fusion protein.

## DRIFT Modeling Parameters
- **Target:** BCR-ABL (constitutively active tyrosine kinase)
- **Mechanism:** ATP-competitive inhibition.
- **DRIFT Config:**
  - `drug_kd`: 10-50 nM (cell-based IC50)
  - `topology`: Modified signaling bridge to include BCR-ABL as the primary driver.

## Predicted Response
1. **Signaling:** Immediate collapse of the downstream RAS/MAPK and PI3K/AKT signaling pathways.
2. **Metabolism:** Shift from high anabolic flux to a quiescent metabolic state.
3. **Phenotype:** High sensitivity in BCR-ABL positive "virtual cells," demonstrating DRIFT's specificity.

## Benchmarking
DRIFT's predictions align with the metabolic reprogramming observed in CML cells upon Imatinib treatment, characterized by decreased glucose uptake and lactate production.
