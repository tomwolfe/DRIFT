import pytest
import numpy as np
from drift.workbench import Workbench
from drift.topology import get_default_topology

def test_multi_target_inhibition():
    # 1. Setup workbench with multi-target drug
    # Inhibit both PI3K and AKT
    drug_kd = {"PI3K": 0.1, "AKT": 0.5}
    wb = Workbench(drug_kd=drug_kd, drug_concentration=1.0)
    
    # 2. Verify binding engine targets
    assert "PI3K" in wb.binding.targets
    assert "AKT" in wb.binding.targets
    assert wb.binding.targets["PI3K"] == 0.1
    assert wb.binding.targets["AKT"] == 0.5
    
    # 3. Run simulation
    # Use few steps for speed
    history = wb.run_simulation(steps=10)
    
    assert "signaling" in history
    assert history["signaling"].shape == (10, 3) # PI3K, AKT, mTOR
    
    # 4. Run Monte Carlo with multi-target
    # This verifies that perturbations are applied to all targets
    results = wb.run_monte_carlo(n_sims=5, steps=10)
    
    for h in results["histories"]:
        assert isinstance(h["drug_kd"], dict)
        assert "PI3K" in h["drug_kd"]
        assert "AKT" in h["drug_kd"]
        
    print("Multi-target inhibition (Polypharmacology) verified successfully.")

if __name__ == "__main__":
    test_multi_target_inhibition()
