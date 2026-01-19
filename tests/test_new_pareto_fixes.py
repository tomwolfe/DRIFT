import pytest
import numpy as np
from drift.workbench import Workbench, clear_worker_cache, _worker_cache
from drift.metabolic import MetabolicBridge, fuzzy_match_reaction_id, DFBASolver

def test_clear_worker_cache():
    # Mock some cache data
    _worker_cache[("test_model", "test_topology", 123)] = "data"
    assert len(_worker_cache) > 0
    
    clear_worker_cache()
    assert len(_worker_cache) == 0

def test_fuzzy_match_reaction_id():
    model_ids = ["EX_glc__D_e", "ATPS4r", "biomass_core"]
    
    # Exact match
    assert fuzzy_match_reaction_id("EX_glc__D_e", model_ids) == "EX_glc__D_e"
    
    # Case insensitive
    assert fuzzy_match_reaction_id("ex_glc__d_e", model_ids) == "EX_glc__D_e"
    
    # Missing prefix
    assert fuzzy_match_reaction_id("glc__D_e", model_ids) == "EX_glc__D_e"
    
    # Double underscore to single (if we implemented it)
    # Our implementation:
    # if "__" in query_id: variations.append(query_id.replace("__", "_"))
    model_ids_single = ["EX_glc_D_e"]
    assert fuzzy_match_reaction_id("EX_glc__D_e", model_ids_single) == "EX_glc_D_e"

def test_bridge_validation_fuzzy():
    class MockReaction:
        def __init__(self, id):
            self.id = id
            
    class MockModel:
        def __init__(self, ids):
            self.reactions = [MockReaction(rid) for rid in ids]

    # Model has EX_glc__D_e
    model = MockModel(["EX_glc__D_e", "BIOMASS_Ecoli_core_w_GAM"])
    
    # Bridge uses glc__D_e (missing EX_)
    bridge = MetabolicBridge(mappings=[
        {"protein_name": "mTOR", "reaction_id": "glc__D_e", "influence": "positive"}
    ])
    
    resolved = bridge.validate_with_model(model)
    assert resolved is True
    assert bridge.mappings[0]["reaction_id"] == "EX_glc__D_e"

def test_workbench_spawn_context_smoke():
    # Just check if it runs without error (smoke test)
    # We use a small model like 'textbook'
    try:
        wb = Workbench(model_name="textbook")
        # If DFBASolver is headless, this might skip some parts, but it should still run
        # We don't need to run a full MC, just check if it initializes
        assert wb.model_name == "textbook"
    except Exception as e:
        pytest.fail(f"Workbench initialization failed: {e}")

if __name__ == "__main__":
    pytest.main([__file__])
