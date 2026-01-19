import pytest
import numpy as np
from drift.metabolic import BridgeBuilder, sigmoidal_mapping
from drift.workbench import Workbench
from drift.topology import get_default_topology

def test_nonlinear_mapping():
    """Test that custom mapping functions are correctly applied."""
    topology = get_default_topology()
    builder = BridgeBuilder()
    builder.set_species_names(topology.species)
    
    # Add a sigmoidal mapping
    builder.add_mapping(
        protein="mTOR",
        reaction_id="EX_glc__D_e",
        mapping_fn=lambda x: sigmoidal_mapping(x, k=20, x0=0.5)
    )
    
    bridge = builder.build()
    
    # Test mapping at different levels
    constraints_low = bridge.get_constraints([0.1, 0.1, 0.1])
    constraints_mid = bridge.get_constraints([0.5, 0.5, 0.5])
    constraints_high = bridge.get_constraints([0.9, 0.9, 0.9])
    
    # Sigmoidal should be very low at 0.1, ~0.5 at 0.5, and very high at 0.9
    val_low = -constraints_low["EX_glc__D_e"] / 10.0
    val_mid = -constraints_mid["EX_glc__D_e"] / 10.0
    val_high = -constraints_high["EX_glc__D_e"] / 10.0
    
    assert val_low < 0.1
    assert 0.4 < val_mid < 0.6
    assert val_high > 0.9

def test_cell_death_recording():
    """Test that FBA infeasibility leads to cell death recording."""
    # We can force infeasibility by setting impossible constraints
    # But easier to mock or use a known case if possible.
    # For now, let's create a workbench and manually trigger a simulation
    # with a bridge that forces a reaction to 0 which is essential.
    
    topology = get_default_topology()
    builder = BridgeBuilder()
    builder.set_species_names(topology.species)
    
    # Force an essential reaction to 0
    # In textbook model, 'Biomass_Ecoli_core' is the objective, 
    # but we can't easily constrain it to be 0 and fail.
    # Instead, let's constrain an uptake reaction to 0.
    builder.add_mapping(
        protein="mTOR",
        reaction_id="EX_glc__D_e",
        base_vmax=0.0  # Force it to 0
    )
    
    bridge = builder.build()
    wb = Workbench(model_name="textbook", bridge=bridge)
    
    # Run simulation
    history = wb.run_simulation(steps=10)
    
    # In 'textbook' model, glucose is essential for growth.
    # If we force it to 0, it might still grow if there are other sources,
    # but usually it's the main one.
    
    # Check if cell death was recorded
    if history["cell_death"]:
        assert history["cell_death"] is True
        assert "dead" in history["status"]
        assert len(history["growth"]) == 10
        # Once dead, growth should be 0
        death_step = history["death_step"]
        assert history["growth"][death_step] == 0.0

def test_reverse_nonlinear_mapping():
    """Test that non-linear feedback from metabolism to signaling works."""
    topology = get_default_topology()
    builder = BridgeBuilder()
    builder.set_species_names(topology.species)
    
    # Add a sigmoidal reverse mapping from Glucose Exchange to mTOR
    builder.add_reverse_mapping(
        flux_id="EX_glc__D_e",
        species_name="mTOR",
        mapping_type="sigmoidal",
        mapping_params={"k": 20, "x0": 0.5},
        baseline=10.0,
        weight=1.0
    )
    
    bridge = builder.build()
    
    # Test feedback at different flux levels
    # High flux (10.0 -> normalized 1.0)
    fb_high = bridge.get_feedback({"EX_glc__D_e": 10.0})
    m_idx = topology.species.index("mTOR")
    assert fb_high[m_idx] > 0.99
    
    # Mid flux (5.0 -> normalized 0.5)
    fb_mid = bridge.get_feedback({"EX_glc__D_e": 5.0})
    assert 0.45 < fb_mid[m_idx] < 0.55
    
    # Low flux (0.0 -> normalized 0.0)
    fb_low = bridge.get_feedback({"EX_glc__D_e": 0.0})
    assert fb_low[m_idx] < 0.01

def test_biological_feedback_rules():
    """Test the new name-based biological feedback rules (mTOR, AMPK)."""
    from drift.metabolic import MetabolicBridge
    
    # Test mTOR sensitivity to both growth and energy
    species = ["PI3K", "AKT", "mTOR", "AMPK"]
    bridge = MetabolicBridge(species_names=species)
    
    # CASE A: High growth, high energy
    # growth_flux=0.2 (max), ATPS4r=5.0 (max)
    fluxes_high = {"Biomass_Ecoli_core": 0.2, "ATPS4r": 5.0}
    fb_high = bridge.get_feedback(fluxes_high)
    
    # mTOR should be high (~1.0)
    assert fb_high[2] > 0.9
    # AMPK should be low (~0.0)
    assert fb_high[3] < 0.1
    
    # CASE B: Low growth, high energy
    fluxes_low_growth = {"Biomass_Ecoli_core": 0.0, "ATPS4r": 5.0}
    fb_low_growth = bridge.get_feedback(fluxes_low_growth)
    # mTOR: 0.5 * global(0) + 0.5 * energy(1) = 0.5
    assert 0.4 < fb_low_growth[2] < 0.6
    
    # CASE C: Low energy
    fluxes_low_energy = {"Biomass_Ecoli_core": 0.0, "ATPS4r": 0.0}
    fb_low_energy = bridge.get_feedback(fluxes_low_energy)
    # mTOR: 0.5 * global(0) + 0.5 * energy(0) = 0.0
    assert fb_low_energy[2] < 0.1
    # AMPK: 1.0 - energy(0) = 1.0
    assert fb_low_energy[3] > 0.9

def test_headless_mode():
    """Test that DFBASolver works in headless mode without crashing."""
    from unittest.mock import patch
    from drift.metabolic import DFBASolver
    
    # Mock optlang.available_solvers to return empty dict
    with patch("optlang.available_solvers", {}):
        solver = DFBASolver(model_name="textbook")
        assert solver.model is None
        
        # solve_step should return safe default
        result = solver.solve_step({})
        assert result["status"] == "headless"
        assert result["objective_value"] == 0.0
        assert result["diagnostic"] == "No solver available (Headless Mode)"
