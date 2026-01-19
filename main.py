from drift.workbench import Workbench
from drift.visualization import create_dashboard
import os
import warnings
import logging

# Suppress cobra and httpx logging noise
logging.getLogger("cobra").setLevel(logging.WARNING)
logging.getLogger("httpx").setLevel(logging.WARNING)

# Suppress cobra warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)

def main():
    print("==========================================================")
    print("   DRIFT: Multi-Scale Stochastic Research Workbench       ")
    print("==========================================================")
    
    # Parameters
    DRUG_KD = 0.5          # High affinity binding (lower is stronger)
    DRUG_CONC = 2.0       # Concentration of drug in system
    SIM_STEPS = 100       # Time-steps for simulation
    MC_ITERATIONS = 30    # Number of Monte Carlo simulations
    
    print(f"[*] Configuration: Kd={DRUG_KD}, [Drug]={DRUG_CONC}")
    print(f"[*] Target Model: Human Metabolic Proxy (textbook-core)")
    
    # Initialize Workbench
    # Note: load_model('textbook') will fetch/load a core metabolic model
    workbench = Workbench(drug_kd=DRUG_KD, drug_concentration=DRUG_CONC, model_name='textbook')
    
    print(f"[*] Running {MC_ITERATIONS} Monte Carlo simulations...")
    all_histories = workbench.run_monte_carlo(n_sims=MC_ITERATIONS, steps=SIM_STEPS)
    
    # Calculate final mean growth to show in CLI
    final_growths = [h['growth'][-1] for h in all_histories]
    avg_final_growth = sum(final_growths) / len(final_growths)
    inhibition = all_histories[0]['inhibition'] * 100
    
    print(f"[+] Target Inhibition: {inhibition:.2f}%")
    print(f"[+] Mean Final Growth Rate: {avg_final_growth:.4f} h⁻¹")
    
    print("[*] Generating interactive multi-scale dashboard...")
    fig = create_dashboard(all_histories)
    
    output_file = "drift_dashboard.html"
    fig.write_html(output_file)
    
    print(f"\n[SUCCESS] Simulation complete.")
    print(f"Results saved to: {os.path.abspath(output_file)}")
    print("==========================================================")

if __name__ == "__main__":
    main()
