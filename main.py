from drift.workbench import Workbench
from drift.visualization import create_dashboard
from drift.config import SimulationConfig
import os
import warnings
import logging
import argparse
import webbrowser
import numpy as np


# Suppress cobra and httpx logging noise
logging.getLogger("cobra").setLevel(logging.WARNING)
logging.getLogger("httpx").setLevel(logging.WARNING)

# Suppress cobra warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)


def main(config: SimulationConfig = None):
    print("==========================================================")
    print("   DRIFT: Multi-Scale Stochastic Research Workbench       ")
    print("==========================================================")

    # Use provided config or create default
    if config is None:
        config = SimulationConfig()

    print(f"[*] Configuration: Kd={config.drug_kd}, [Drug]={config.drug_concentration}")
    print(f"[*] Target Model: Human Metabolic Proxy ({config.model_name})")

    try:
        # Initialize Workbench
        # Note: load_model('textbook') will fetch/load a core metabolic model
        workbench = Workbench(
            drug_kd=config.drug_kd,
            drug_concentration=config.drug_concentration,
            model_name=config.model_name,
        )

        print(f"[*] Running {config.mc_iterations} Monte Carlo simulations...")
        results = workbench.run_monte_carlo(
            n_sims=config.mc_iterations, steps=config.sim_steps, n_jobs=config.n_jobs
        )
        all_histories = results["histories"]
        basal_growth = results["basal_growth"]

        print(f"[+] Completed {len(all_histories)} simulations")

        # Calculate final mean growth to show in CLI
        final_growths = [h["growth"][-1] for h in all_histories]
        avg_final_growth = sum(final_growths) / len(final_growths)
        std_final_growth = np.std(final_growths)
        inhibition = all_histories[0]["inhibition"] * 100

        # Calculate relative growth (normalization)
        rel_growth = (avg_final_growth / basal_growth) * 100 if basal_growth > 0 else 0

        print("\n--- Simulation Summary ---")
        print(f"Target Inhibition:      {inhibition:.2f}%")
        print(f"Basal Growth Rate:      {basal_growth:.4f} h⁻¹")
        print(
            f"Mean Treated Growth:    {avg_final_growth:.4f} h⁻¹ (±{std_final_growth:.4f})"
        )
        print(f"Relative Vitality:      {rel_growth:.1f}%")
        print("--------------------------\n")

        print("[*] Generating interactive multi-scale dashboard...")
        fig = create_dashboard(results)

        output_file = "drift_dashboard.html"
        fig.write_html(output_file)

        print("[SUCCESS] Simulation complete.")
        print(f"Results saved to: {os.path.abspath(output_file)}")

        # Pareto Optimal: Auto-open the dashboard to reduce user friction
        print("[*] Opening dashboard in browser...")
        webbrowser.open(f"file://{os.path.abspath(output_file)}")
        print("==========================================================")

    except Exception as e:
        print(f"[ERROR] Simulation failed: {str(e)}")
        raise


def parse_args():
    parser = argparse.ArgumentParser(
        description="DRIFT: Multi-Scale Stochastic Research Workbench"
    )
    parser.add_argument("--config", type=str, help="Path to JSON configuration file")
    parser.add_argument("--drug-kd", type=float, help="Drug dissociation constant")
    parser.add_argument("--drug-conc", type=float, help="Drug concentration")
    parser.add_argument("--sim-steps", type=int, help="Number of simulation steps")
    parser.add_argument(
        "--mc-iterations", type=int, help="Number of Monte Carlo iterations"
    )
    parser.add_argument("--model-name", type=str, help="Metabolic model name")
    parser.add_argument(
        "--output", type=str, default="drift_dashboard.html", help="Output file name"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Load configuration from file or command line arguments
    if args.config:
        config = SimulationConfig.from_json_file(args.config)
    else:
        # Create config with default values, overriding with any provided args
        config_kwargs = {}
        if args.drug_kd is not None:
            config_kwargs["drug_kd"] = args.drug_kd
        if args.drug_conc is not None:
            config_kwargs["drug_concentration"] = args.drug_conc
        if args.sim_steps is not None:
            config_kwargs["sim_steps"] = args.sim_steps
        if args.mc_iterations is not None:
            config_kwargs["mc_iterations"] = args.mc_iterations
        if args.model_name is not None:
            config_kwargs["model_name"] = args.model_name

        config = SimulationConfig(**config_kwargs)

    main(config)
