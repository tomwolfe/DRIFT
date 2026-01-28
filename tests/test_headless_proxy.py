import pytest
import numpy as np
from drift.metabolic import DFBASolver, MetabolicBridge

def test_limiting_substrate_proxy():
    # 1. Setup headless solver
    solver = DFBASolver(model_name="textbook")
    solver.headless = True
    solver.model = None # Force proxy logic
    solver.proxy_params = {
        "base_growth": 0.2,
        "yield_coefficient": 0.02
    }
    assert solver.headless == True
    
    # 2. Setup bridge with two mappings
    # One is highly restrictive (low scaling), one is not.
    bridge = MetabolicBridge(
        species_names=["A", "B"],
        mappings=[
            {
                "protein_name": "A",
                "reaction_id": "REAC1",
                "influence": "positive",
                "base_vmax": 10.0,
                "basal_scaling": 0.1,
                "max_scaling": 0.9
            },
            {
                "protein_name": "B",
                "reaction_id": "REAC2",
                "influence": "positive",
                "base_vmax": 10.0,
                "basal_scaling": 0.1,
                "max_scaling": 0.9
            }
        ]
    )
    
    # Case A: A is very low (bottleneck)
    # A=0.0 -> scaling_A = 0.1
    # B=1.0 -> scaling_B = 1.0
    state_a = np.array([0.0, 1.0])
    constraints_a, scalings_a = bridge.get_constraints_with_scalings(state_a)
    result_a = solver.solve_step(constraints_a, scalings=scalings_a)
    growth_a = result_a["objective_value"]
    
    # Case B: Both are high
    # A=1.0 -> scaling_A = 1.0
    # B=1.0 -> scaling_B = 1.0
    state_b = np.array([1.0, 1.0])
    constraints_b, scalings_b = bridge.get_constraints_with_scalings(state_b)
    result_b = solver.solve_step(constraints_b, scalings=scalings_b)
    growth_b = result_b["objective_value"]
    
    # Case C: B is very low (bottleneck)
    # A=1.0 -> scaling_A = 1.0
    # B=0.0 -> scaling_B = 0.1
    state_c = np.array([1.0, 0.0])
    constraints_c, scalings_c = bridge.get_constraints_with_scalings(state_c)
    result_c = solver.solve_step(constraints_c, scalings=scalings_c)
    growth_c = result_c["objective_value"]
    
    # Growth should be much lower in case A and C than in B
    assert growth_a < growth_b
    assert growth_c < growth_b
    
    # In this simple model, growth_a and growth_c should be similar because they have same min_scaling
    assert np.isclose(growth_a, growth_c)
    
    # Check if growth is proportional to scaling^0.7 as implemented
    # scaling = 0.1 -> growth = 0.2 * (0.1**0.7) = 0.2 * 0.199 = 0.0399
    # scaling = 1.0 -> growth = 0.2 * (1.0**0.7) = 0.2
    assert np.isclose(growth_b, 0.2)
    assert growth_a < 0.05
    
    print("Limiting substrate proxy verified successfully.")

if __name__ == "__main__":
    test_limiting_substrate_proxy()
