# 🖼️ Case Study Gallery

The DRIFT Gallery showcases the framework's ability to predict multi-scale responses for a wide range of clinically relevant compounds. Each case study integrates known biochemical data with DRIFT's stochastic signaling and metabolic engines.

## 💊 Featured Case Studies

| Drug | Target | Primary Indication | Case Study Link |
| --- | --- | --- | --- |
| **Alpelisib** | PI3Kα | Breast Cancer | [View Case Study](validation.md#9-experimental-benchmarking-alpelisib-byl719-case-study) |
| **Erlotinib** | EGFR | NSCLC | [View Case Study](gallery/erlotinib.md) |
| **Imatinib** | BCR-ABL | CML | [View Case Study](gallery/imatinib.md) |
| **Everolimus** | mTOR | Various Cancers | [View Case Study](gallery/everolimus.md) |
| **Cetuximab** | EGFR (mAb) | Colorectal Cancer | [View Case Study](gallery/cetuximab.md) |
| **Palbociclib** | CDK4/6 | Breast Cancer | [View Case Study](gallery/palbociclib.md) |

## 🧪 About the Gallery

These case studies are generated using the DRIFT `Workbench`. They demonstrate how:
1.  **Binding Affinity ($K_d$):** Translates to signaling inhibition.
2.  **Signaling Dynamics:** Captures the temporal response of key nodes like AKT and mTOR.
3.  **Metabolic Phenotype:** Predicts growth inhibition and flux redistribution.

### 🚀 Running Your Own Comparison
You can reproduce these results or test new drugs using the `examples/drug_comparison.py` script.
