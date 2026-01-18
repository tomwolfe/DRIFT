# System Architecture

DRIFT employs a hierarchical coupling strategy to link three distinct biological time-scales and resolutions.

## 1. The Molecular Scale (Seconds)
**Module:** `BindingEngine`
- **Logic:** Equilibrium binding kinetics.
- **Implementation:** The Hill Equation determines the fraction of target proteins bound by the ligand. This acts as a continuous "throttle" on the signaling cascade.
- **Formula:** $\theta = \frac{[D]^n}{[D]^n + K_d^n}$

## 2. The Signaling Scale (Minutes)
**Module:** `StochasticIntegrator`
- **Logic:** Langevin Dynamics (SDEs).
- **Implementation:** We model the PI3K/AKT/mTOR axis as a series of coupled differential equations with an additive Gaussian noise term $\eta(t)$.
- **Coupling:** The `inhibition` parameter from the Molecular Scale directly reduces the rate of PI3K-mediated AKT activation.
- **Stochasticity:** Captures the "intrinsic noise" of the cytoplasm, ensuring that no two cells (simulations) follow the exact same trajectory.

## 3. The Metabolic Scale (Hours)
**Module:** `DFBASolver` & `MetabolicBridge`
- **Logic:** Constraint-Based Modeling (CBM).
- **Implementation:** Flux Balance Analysis (FBA) assumes a pseudo-steady state for internal metabolites while maximizing an objective function (e.g., Biomass production).
- **Coupling:** The `MetabolicBridge` maps the concentration of `mTOR` from the Signaling Scale to the `Vmax` (upper/lower bounds) of specific metabolic reactions, such as Glucose uptake.

## Multi-Scale Integration Loop
1. **Binding:** Calculate fixed inhibition level.
2. **Signaling Step:** Update protein states $S_{t} \to S_{t+dt}$ using Langevin integration.
3. **Bridge:** Translate $S_{t+dt}$ to metabolic bounds $B_{t+dt}$.
4. **Metabolic Step:** Solve $max(f \cdot v)$ subject to $S \cdot v = 0$ and $v < B_{t+dt}$.
5. **Record:** Store state and flux distribution.
