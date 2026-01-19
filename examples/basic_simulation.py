from drift.workbench import Workbench
import matplotlib.pyplot as plt
import numpy as np

def main():
    """
    This example demonstrates how to use the DRIFT Workbench to run a single
    multi-scale simulation of drug-induced metabolic drift.
    """
    print("Initializing DRIFT Workbench...")

    # Parameters
    DRUG_KD = 0.5
    DRUG_CONCENTRATION = 1.0
    STEPS = 200

    # 1. Initialize the Workbench
    # This loads the metabolic model and prepares the engines
    wb = Workbench(drug_kd=DRUG_KD, drug_concentration=DRUG_CONCENTRATION, model_name='textbook')

    # 2. Run a single temporal simulation
    # This performs the integrated Binding -> Signaling -> Metabolic loop
    print(f"Running simulation for {STEPS} steps...")
    history = wb.run_simulation(steps=STEPS)

    # 3. Analyze and visualize the results
    time = history['time']
    signaling = history['signaling']  # PI3K, AKT, mTOR
    growth = history['growth']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Plot Signaling Trajectories
    ax1.plot(time, signaling[:, 0], label='PI3K (Inhibited)')
    ax1.plot(time, signaling[:, 1], label='AKT')
    ax1.plot(time, signaling[:, 2], label='mTOR')
    ax1.set_ylabel('Relative Concentration')
    ax1.set_title(f'Signaling Dynamics (Drug Conc: {DRUG_CONCENTRATION})')
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.7)

    # Plot Metabolic Response
    ax2.plot(time, growth, color='forestgreen', linewidth=2)
    ax2.set_ylabel('Growth Rate (hr^-1)')
    ax2.set_xlabel('Time (arbitrary units)')
    ax2.set_title('Metabolic Output (Growth Rate)')
    ax2.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    print("Simulation complete. Displaying results (if a display is available)...")
    # plt.show() # Uncomment to see the plot

    # Print summary statistics
    print("\nSimulation Summary:")
    print(f"Initial Growth: {growth[0]:.4f}")
    print(f"Final Growth (Mean of last 10): {np.mean(growth[-10:]):.4f}")
    print(f"Binding Inhibition: {history['inhibition']*100:.1f}%")

if __name__ == "__main__":
    main()
