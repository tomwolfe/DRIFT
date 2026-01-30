import numpy as np
import pytest
from drift.workbench import Workbench
from drift.metabolic import MetabolicBridge, BridgeBuilder
from drift.topology import get_default_topology
from drift.core.exceptions import DesyncError

def test_reproducibility_parallel_seeds():
    """
    Spawns 4 Workbench instances with the same seed and asserts that
    signaling trajectories are bit-perfectly identical.
    """
    seed = 42
    steps = 50
    n_instances = 4
    
    histories = []
    for i in range(n_instances):
        wb = Workbench(random_state=seed)
        hist = wb.run_simulation(steps=steps)
        histories.append(hist)
        
    # Assert all signaling trajectories are identical
    base_signaling = histories[0]["signaling"]
    for i in range(1, n_instances):
        assert np.array_equal(base_signaling, histories[i]["signaling"]), f"Instance {i} signaling trajectory differs from instance 0"
    
    # Assert growth curves are also identical
    base_growth = histories[0]["growth"]
    for i in range(1, n_instances):
        assert np.array_equal(base_growth, histories[i]["growth"]), f"Instance {i} growth curve differs from instance 0"

def test_desync_validation():
    """
    Verifies that Workbench raises DesyncError if the bridge mappings 
    refer to species not present in the topology.
    """
    topology = get_default_topology()
    # Create a bridge with a non-existent protein
    bridge = BridgeBuilder().add_mapping(protein_name="NON_EXISTENT_PROTEIN", reaction_id="EX_glc__D_e").build()
    
    with pytest.raises(DesyncError) as excinfo:
        Workbench(topology=topology, bridge=bridge)
    
    assert "Structural Integrity Failure" in str(excinfo.value)
    assert "NON_EXISTENT_PROTEIN" in str(excinfo.value)

def test_proxy_hardening_bounds():
    """
    Harden the proxy by ensuring it raises SolverError on physical bound violations.
    """
    from drift.metabolic import DFBASolver
    from drift.core.exceptions import SolverError
    solver = DFBASolver(model_name="textbook")
    # Manually force headless state for this instance to test proxy logic
    solver.headless = True
    if not hasattr(solver, "proxy_params"):
        solver.proxy_params = {"base_growth": 0.2}

    # Test growth > 2.0 hr⁻¹ (should raise SolverError)
    solver.proxy_params["base_growth"] = 3.0
    with pytest.raises(SolverError):
        solver._solve_proxy(constraints={"EX_glc__D_e": -10.0})

    # Test growth < 0.0
    solver.proxy_params["base_growth"] = -0.1
    with pytest.raises(SolverError):
        solver._solve_proxy(constraints={"EX_glc__D_e": -10.0})

def test_local_rng_independence():
    """
    Ensures that two Workbench instances with different seeds remain independent.
    """
    wb1 = Workbench(random_state=1)
    wb2 = Workbench(random_state=2)
    
    h1 = wb1.run_simulation(steps=20)
    h2 = wb2.run_simulation(steps=20)
    
    assert not np.array_equal(h1["signaling"], h2["signaling"])
