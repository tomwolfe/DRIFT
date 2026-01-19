import pytest
import numpy as np
from drift.workbench import Workbench
from drift.topology import Topology
from drift.metabolic import BridgeBuilder

def test_custom_topology_and_bridge_builder():
    """Test that a user can define a custom topology and use BridgeBuilder."""
    # 1. Define custom topology (MAPK-like)
    custom_species = ["RAS", "RAF", "MEK", "ERK"]
    custom_params = {"k_deg": 0.1}
    
    # Simple drift: each activates the next
    def custom_drift(state, params, inhibition):
        ras, raf, mek, erk = state
        # Drug inhibits RAS
        effective_ras = ras * (1.0 - inhibition)
        
        dras = 0.1 - 0.1 * ras
        draf = 0.5 * effective_ras * (1.0 - raf) - 0.1 * raf
        dmek = 0.5 * raf * (1.0 - mek) - 0.1 * mek
        derk = 0.5 * mek * (1.0 - erk) - 0.1 * erk
        
        return np.array([dras, draf, dmek, derk])

    topology = Topology(
        species=custom_species,
        parameters=custom_params,
        drift_fn=custom_drift,
        name="MAPK"
    )

    # 2. Use BridgeBuilder to link ERK (idx 3) to glucose uptake
    bridge = (
        BridgeBuilder()
        .add_mapping(protein_idx=3, reaction_id="EX_glc__D_e", influence="positive", base_vmax=15.0)
        .build()
    )

    # 3. Initialize Workbench with custom components
    wb = Workbench(
        drug_kd=0.1,
        drug_concentration=1.0,
        model_name="textbook",
        topology=topology,
        bridge=bridge
    )

    # 4. Run simulation
    history = wb.run_simulation(steps=10)
    
    assert "signaling" in history
    assert history["signaling"].shape == (10, 4)
    assert "growth" in history
    assert len(history["growth"]) == 10
    assert history["inhibition"] > 0

def test_bridge_builder_validation():
    """Test BridgeBuilder validation."""
    builder = BridgeBuilder()
    with pytest.raises(ValueError, match="influence must be 'positive' or 'negative'"):
        builder.add_mapping(0, "RXN", influence="maybe")

def test_workbench_invalid_inputs():
    """Test Workbench raises errors for invalid inputs."""
    with pytest.raises(ValueError, match="drug_kd must be positive"):
        Workbench(drug_kd=-1)
    
    wb = Workbench()
    with pytest.raises(ValueError, match="steps must be positive"):
        wb.run_simulation(steps=0)
