import pytest
import numpy as np
from drift.topology import Topology, drift_model
from drift.binding import BindingEngine
from drift.metabolic import MetabolicBridge, DFBASolver
from drift.engine import SimulationEngine
from drift.visualization import create_dashboard

@drift_model("Synthetic_Circuit")
def synthetic_drift(state, params, feedback=None):
    # 2-node synthetic circuit: A -> B
    # A is produced at constant rate, B is produced by A
    a, b = state[0], state[1]
    
    # params: [k_a_prod, k_a_deg, k_b_prod, k_b_deg, inh_a, inh_b]
    ka_p, ka_d = params[0], params[1]
    kb_p, kb_d = params[2], params[3]
    
    # Handle multi-target inhibition
    if len(params) == 5:
        inh_a = params[4]
        inh_b = 0.0
    elif len(params) >= 6:
        inh_a, inh_b = params[4], params[5]
    else:
        inh_a = inh_b = 0.0
        
    fb = feedback if feedback is not None else np.ones(2)
    
    da = ka_p * fb[0] * (1.0 - inh_a) - ka_d * a
    db = kb_p * fb[1] * a * (1.0 - inh_b) - kb_d * b
    
    res = np.zeros_like(state)
    res[0], res[1] = da, db
    return res

def test_generalized_topology_dashboard():
    # 1. Setup generalized 2-node topology
    topo = Topology(
        species=["NodeA", "NodeB"],
        parameters={
            "ka_p": 0.2, "ka_d": 0.1,
            "kb_p": 0.5, "kb_d": 0.1
        },
        drift_fn=synthetic_drift,
        name="Synthetic_2Node"
    )
    
    from drift.signaling import StochasticIntegrator
    integrator = StochasticIntegrator(dt=0.1, noise_scale=0.01, topology=topo)
    
    # 2. Setup headless solver and bridge
    solver = DFBASolver(model_name="textbook") # Will be headless if no solver
    bridge = MetabolicBridge(
        species_names=["NodeA", "NodeB"],
        mappings=[
            {
                "protein_name": "NodeB",
                "reaction_id": "EX_glc__D_e",
                "influence": "positive",
                "base_vmax": 10.0
            }
        ]
    )
    
    # 3. Setup binding with multi-target inhibition
    # Inhibit both NodeA and NodeB
    binding = BindingEngine(targets={"NodeA": 0.5, "NodeB": 1.5})
    
    engine = SimulationEngine(
        integrator=integrator,
        solver=solver,
        bridge=bridge,
        binding=binding
    )
    
    # 4. Run simulation
    result = engine.run(steps=50, drug_concentration=1.0)
    
    assert "species_names" in result
    assert result["species_names"] == ["NodeA", "NodeB"]
    assert result["signaling"].shape == (50, 2)
    
    # 5. Verify dashboard generation
    # create_dashboard expects a list of histories
    dashboard_input = {"histories": [result], "basal_growth": 0.2}
    fig = create_dashboard(dashboard_input)
    
    assert fig is not None
    # Check if the trace names match NodeA and NodeB
    trace_names = [t.name for t in fig.data if t.name in ["NodeA", "NodeB"]]
    assert "NodeA" in trace_names
    assert "NodeB" in trace_names
    
    print("Generalized topology and dashboard verified successfully.")

if __name__ == "__main__":
    test_generalized_topology_dashboard()
