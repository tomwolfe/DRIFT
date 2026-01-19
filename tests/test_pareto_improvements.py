
import numpy as np
import pytest
from drift.workbench import Workbench
from drift.topology import get_complex_topology
from drift.metabolic import BridgeBuilder

def test_reproducibility():
    """Verify that random_state and per-simulation seeds ensure bit-perfect reproducibility."""
    topo = get_complex_topology()
    wb = Workbench(topology=topo, random_state=42)
    
    # Run two simulations with same seed
    res1 = wb.run_simulation(steps=20, seed=123)
    res2 = wb.run_simulation(steps=20, seed=123)
    
    # Run one with different seed
    res3 = wb.run_simulation(steps=20, seed=456)
    
    assert np.allclose(res1["signaling"], res2["signaling"]), "Identical seeds produced different results"
    assert not np.allclose(res1["signaling"], res3["signaling"]), "Different seeds produced identical results"

def test_complex_feedback_probes():
    """Verify that flux probes correctly drive signaling nodes (Mechanistic Feedback)."""
    topo = get_complex_topology()
    
    # Build a bridge that uses the 'ATP_total' probe for AMPK (index 3)
    bridge = BridgeBuilder() \
        .set_species_names(topo.species) \
        .add_reverse_mapping(flux_id="ATP_total", species_name="AMPK", influence="negative", weight=0.8, baseline=25.0) \
        .build()
        
    wb = Workbench(topology=topo, bridge=bridge, random_state=42)
    res = wb.run_simulation(steps=10)
    
    assert "signaling" in res
    assert res["signaling"].shape[1] == 4, "Topology should have 4 species"
    assert "growth" in res
    assert len(res["growth"]) == 10

def test_schema_enforcement():
    """Verify that simulation results follow the SimulationResult TypedDict schema."""
    wb = Workbench(random_state=42)
    res = wb.run_simulation(steps=5)
    
    # Check required fields from TypedDict
    required_keys = [
        "time", "signaling", "growth", "status", "cell_death", 
        "death_step", "death_cause", "inhibition", "drug_kd", 
        "params", "headless"
    ]
    for key in required_keys:
        assert key in res, f"Missing required key in simulation result: {key}"
    
    assert isinstance(res["time"], np.ndarray)
    assert isinstance(res["signaling"], np.ndarray)
    assert isinstance(res["growth"], np.ndarray)

def test_solver_aware_proxy(monkeypatch):
    """Verify that headless mode uses the enhanced proxy logic."""
    from drift.metabolic import DFBASolver
    
    # Mock _check_solver to return empty list, forcing headless mode
    monkeypatch.setattr(DFBASolver, "_check_solver", lambda self: [])
    
    solver = DFBASolver(model_name="non_existent", strict=False)
    assert solver.headless
    
    constraints = {"EX_glc__D_e": -10.0}
    result = solver.solve_step(constraints)
    
    assert result["status"] == "proxied"
    assert "QUALITATIVE_NOTICE" in result
    assert "Solver-Aware Proxy" in result["diagnostic"]
    # Check that energy proxy is present
    assert "ATPS4r" in result["fluxes"]
