"""
Example: Custom Signaling Pathway and Metabolic Bridge

This example demonstrates how to:
1. Define a custom signaling topology (MAPK pathway).
2. Create a custom metabolic bridge using the BridgeBuilder.
3. Run a multi-scale simulation linking the two.
"""

import numpy as np
import matplotlib.pyplot as plt
from drift.workbench import Workbench
from drift.topology import Topology
from drift.metabolic import BridgeBuilder

def main():
    print("Initializing custom MAPK signaling topology...")
    
    # 1. Define custom topology species
    species = ["RAS", "RAF", "MEK", "ERK"]
    
    # 2. Define kinetic parameters
    parameters = {
        "k_basal": 0.1,
        "k_act": 0.8,
        "k_deact": 0.2,
    }

    # 3. Define drift function (mechanistic SDE components)
    def mapk_drift(state, params, inhibition):
        ras, raf, mek, erk = state
        
        # Drug inhibits RAS activation
        effective_ras = ras * (1.0 - inhibition)
        
        # RAS basal dynamics
        dras = 0.1 - 0.1 * ras
        # Cascade: RAS -> RAF -> MEK -> ERK
        draf = 0.8 * effective_ras * (1.0 - raf) - 0.2 * raf
        dmek = 0.8 * raf * (1.0 - mek) - 0.2 * mek
        derk = 0.8 * mek * (1.0 - erk) - 0.2 * erk
        
        return np.array([dras, draf, dmek, derk])

    topology = Topology(
        species=species,
        parameters=parameters,
        drift_fn=mapk_drift,
        name="MAPK_Cascade"
    )

    print("Building custom metabolic bridge...")
    # 4. Use BridgeBuilder to link ERK activity to glucose uptake
    # ERK (index 3) is a known activator of GLUT1 translocation and hexokinase activity
    bridge = (
        BridgeBuilder()
        .add_mapping(
            protein_idx=3, 
            reaction_id="EX_glc__D_e", 
            influence="positive", 
            base_vmax=12.0
        )
        .build()
    )

    print("Running Workbench simulation...")
    # 5. Initialize Workbench
    wb = Workbench(
        drug_kd=0.05, 
        drug_concentration=0.5, 
        model_name="textbook", 
        topology=topology,
        bridge=bridge
    )

    # 6. Run simulation
    steps = 150
    history = wb.run_simulation(steps=steps)

    # 7. Visualize results
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    time = history["time"]
    
    # Plot signaling
    for i, name in enumerate(species):
        ax1.plot(time, history["signaling"][:, i], label=name)
    ax1.set_ylabel("Normalized Concentration")
    ax1.set_title(f"Signaling Dynamics (Inhibition: {history['inhibition']:.2f})")
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Plot growth
    ax2.plot(time, history["growth"], color="black", linewidth=2)
    ax2.set_ylabel("Growth Rate (hr^-1)")
    ax2.set_xlabel("Time (arbitrary units)")
    ax2.set_title("Metabolic Phenotype (Growth)")
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    print("Simulation complete. Saving plot to 'custom_pathway_output.png'...")
    plt.savefig("outputs/custom_pathway_output.png")
    print("Done!")

if __name__ == "__main__":
    main()
