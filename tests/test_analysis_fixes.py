import pytest
import numpy as np
import os
import pandas as pd
from numba import njit
from drift.signaling import langevin_step, create_langevin_integrator
from drift.metabolic import MetabolicBridge, BridgeBuilder, DFBASolver
from drift.workbench import Workbench
from drift.validation import ExternalValidator
from drift.core.exceptions import DesyncError

@njit
def extreme_negative_drift(state, params, feedback=None):
    """A drift function that always pushes the state strongly negative."""
    # Numba needs this to be a top-level jitted function for create_langevin_integrator
    return -100.0 * np.ones_like(state)

def test_sde_clamping():
    """Test that SDE state is clamped to [0, 1]."""
    # We'll use a very large noise and dt to force boundary crossing
    state = np.array([0.05, 0.5, 0.95])
    dt = 1.0
    params = np.array([0.1, 0.1, 0.5, 0.1, 0.5, 0.1, 0.0]) # PI3K_AKT_mTOR params
    noise_scale = 5.0 # Huge noise
    
    for _ in range(100):
        state = langevin_step(state, dt, params, noise_scale)
        assert np.all(state >= 0)
        assert np.all(state <= 1)
        # Ensure it's exactly 0 or 1 sometimes with this much noise
        if np.any((state == 0) | (state == 1)):
            break
    else:
        # Should have hit a boundary with noise_scale=5.0
        pass

def test_strict_mapping():
    """Test that strict mapping raises errors/warnings correctly."""
    builder = BridgeBuilder()
    builder.add_mapping(protein="mTOR", reaction_id="NON_EXISTENT_RXN")
    builder.set_strict_mapping(True)
    bridge = builder.build()
    
    solver = DFBASolver(model_name="textbook")
    if not solver.headless:
        # Should return False because NON_EXISTENT_RXN isn't there
        assert bridge.validate_with_model(solver.model) is False

def test_configurable_scaling():
    """Test that scaling constants are respected."""
    builder = BridgeBuilder()
    builder.set_species_names(["PI3K", "AKT", "mTOR"])
    # Basal 0.5, Max 0.5 -> Range [0.5, 1.0]
    builder.add_mapping(protein="mTOR", reaction_id="EX_glc__D_e", basal_scaling=0.5, max_scaling=0.5)
    bridge = builder.build()
    
    # If mTOR is 0, scaling should be 0.5
    constraints = bridge.get_constraints(np.array([0.0, 0.0, 0.0]))
    assert constraints["EX_glc__D_e"] == -5.0 # - (10.0 * 0.5)
    
    # If mTOR is 1.0, scaling should be 1.0
    constraints = bridge.get_constraints(np.array([0.0, 0.0, 1.0]))
    assert constraints["EX_glc__D_e"] == -10.0 # - (10.0 * 1.0)

def test_auto_generator_mode():
    """Test that Workbench automatically switches to generator for large runs."""
    workbench = Workbench(model_name="textbook")
    gen = workbench.run_monte_carlo(n_sims=2, steps=10, return_generator=True)
    assert hasattr(gen, "__iter__")
    results = list(gen)
    assert len(results) == 2

def test_external_validator(tmp_path):
    """Test the external validation framework."""
    # Create dummy experimental data
    df = pd.DataFrame({
        "time": [0, 5, 10],
        "AKT": [0.1, 0.5, 0.9],
        "growth": [0.2, 0.15, 0.1]
    })
    csv_path = tmp_path / "exp_data.csv"
    df.to_csv(csv_path, index=False)
    
    validator = ExternalValidator(str(csv_path))
    
    # Dummy simulation history
    hist = {
        "time": np.array([0, 2, 4, 6, 8, 10]),
        "signaling": np.zeros((6, 3)),
        "growth": np.array([0.2, 0.19, 0.18, 0.17, 0.16, 0.15])
    }
    hist["signaling"][:, 1] = np.linspace(0.1, 0.9, 6) # AKT
    
    mapping = {"AKT": "AKT", "growth": "growth"}
    metrics = validator.benchmark(hist, mapping)
    
    assert "AKT" in metrics
    assert "growth" in metrics
    assert metrics["AKT"]["correlation"] > 0.99

def test_desync_error_missing_species_names():
    """Test DesyncError when species_names are missing."""
    bridge = MetabolicBridge(mappings=[{"protein_name": "mTOR", "reaction_id": "EX_glc__D_e"}])
    # bridge.species_names will be None because mappings is provided but species_names is not
    with pytest.raises(DesyncError, match="requires species_names"):
        bridge.get_constraints(np.array([0.5]))

def test_desync_error_length_mismatch():
    """Test DesyncError when state length doesn't match species_names."""
    bridge = MetabolicBridge(species_names=["PI3K", "AKT", "mTOR"])
    with pytest.raises(DesyncError, match="desynchronized"):
        bridge.get_constraints(np.array([0.5, 0.5])) # Length 2, expected 3

def test_desync_error_missing_protein_in_state():
    """Test DesyncError when a protein in mapping is not in species_names."""
    builder = BridgeBuilder()
    builder.set_species_names(["PI3K", "AKT"])
    builder.add_mapping(protein="mTOR", reaction_id="EX_glc__D_e")
    bridge = builder.build()
    with pytest.raises(DesyncError, match="not found in signaling state"):
        bridge.get_constraints(np.array([0.5, 0.5]))

def test_strict_type_checking():
    """Test that get_constraints enforces numpy array or list/tuple."""
    bridge = MetabolicBridge(species_names=["mTOR"])
    with pytest.raises(ValueError, match="must be a numpy array"):
        bridge.get_constraints("not an array")
    
    # List should be converted and work
    constraints = bridge.get_constraints([0.5])
    assert "EX_glc__D_e" in constraints

def test_dfba_error_fingerprint():
    """Test that DFBASolver returns an error_fingerprint on infeasibility."""
    solver = DFBASolver(model_name="textbook")
    if solver.headless:
        pytest.skip("Test requires a real LP solver")
    
    # Force infeasibility by setting impossible bounds
    result = solver.optimize({"EX_glc__D_e": 10.0}) # Positive LB for uptake is usually infeasible
    assert result["status"] != "optimal"
    assert "error_fingerprint" in result
    assert len(result["error_fingerprint"]) == 8

def test_dfba_solve_step_alias():
    """Test that solve_step still works as an alias for optimize."""
    solver = DFBASolver(model_name="textbook")
    constraints = {"EX_glc__D_e": -5.0}
    result = solver.solve_step(constraints)
    assert "status" in result

def test_milstein_non_negativity_clamping():
    """Test the 'if val < 0: val = 0' logic specifically."""
    step_fn = create_langevin_integrator(extreme_negative_drift)
    
    state = np.array([0.1, 0.1, 0.1])
    dt = 0.1
    # params length must match what step_fn expects (len(species) + len(species))
    # For PI3K/AKT/mTOR default it's 3 + 3 = 6 or similar.
    # We'll just provide a large enough array.
    params = np.zeros(20)
    noise_scale = 0.0
    
    new_state = step_fn(state, dt, params, noise_scale)
    assert np.all(new_state >= 0)
    assert np.all(new_state <= 1)
    # With -100 drift and 0.1 dt, 0.1 - 10.0 = -9.9, should be clamped to 0.
    assert np.all(new_state == 0)

def test_desync_error_invalid_mapping_config():
    """Test DesyncError when mapping has neither name nor index."""
    bridge = MetabolicBridge(species_names=["mTOR"])
    # Manually corrupt mappings
    bridge.mappings = [{"reaction_id": "EX_glc__D_e"}] 
    with pytest.raises(DesyncError, match="lacks both valid protein_name and protein_idx"):
        bridge.get_constraints(np.array([0.5]))

def test_metabolic_bridge_calibration():
    """Test bridge calibration logic."""
    bridge = MetabolicBridge()
    original_basal = bridge.basal_growth_rate
    bridge.calibrate({"Biomass_Ecoli_core": 0.5})
    assert bridge.basal_growth_rate == 0.5
    assert bridge.basal_growth_rate != original_basal

def test_michaelis_menten_mapping():
    """Test Michaelis-Menten mapping function."""
    from drift.metabolic import michaelis_menten_mapping
    # At x=1.0, should return 1.0 due to normalization in our implementation
    assert michaelis_menten_mapping(1.0, km=0.5) == 1.0
    assert michaelis_menten_mapping(0.0, km=0.5) == 0.0
    # At x=km, normalized MM is (km*(km+1))/(km+km) = (km+1)/2
    assert michaelis_menten_mapping(0.5, km=0.5) == 0.75

def test_sigmoidal_mapping():
    """Test Sigmoidal mapping function."""
    from drift.metabolic import sigmoidal_mapping
    # At x=x0, should return 0.5
    assert sigmoidal_mapping(0.5, k=10, x0=0.5) == 0.5
    assert sigmoidal_mapping(1.0, k=100, x0=0.5) > 0.99
    assert sigmoidal_mapping(0.0, k=100, x0=0.5) < 0.01

def test_fuzzy_match_reaction_id():
    """Test fuzzy matching for reaction IDs."""
    from drift.metabolic import fuzzy_match_reaction_id
    model_ids = ["EX_glc__D_e", "ATPS4r", "PGI"]
    # Exact and prefix matches should work
    assert fuzzy_match_reaction_id("glc__D_e", model_ids) == "EX_glc__D_e"
    assert fuzzy_match_reaction_id("EX_glc__D_e", model_ids) == "EX_glc__D_e"
    # Case insensitive
    assert fuzzy_match_reaction_id("EX_GLC__D_E", model_ids) == "EX_glc__D_e"
    # Similarity
    assert fuzzy_match_reaction_id("ATPS4r", model_ids) == "ATPS4r"
    assert fuzzy_match_reaction_id("NOT_THERE", model_ids) is None
