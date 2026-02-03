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


def main(config: SimulationConfig = None, run_sensitivity=False, json_export=None, human_cancer=False, diabetes=False, web_schema=None):
    print("==========================================================")
    print("   DRIFT: Multi-Scale Stochastic Research Workbench       ")
    print("==========================================================")

    # Use provided config or create default
    if config is None:
        config = SimulationConfig()
    
    # Export web schema if requested
    if web_schema:
        import json
        schema = config.to_web_schema()
        with open(web_schema, "w") as f:
            json.dump(schema, f, indent=4)
        print(f"[*] Web configuration schema exported to: {web_schema}")
        return

    # Apply presets if requested
    topology = None
    bridge = None
    if human_cancer:
        from drift.presets import get_human_cancer_topology, get_human_cancer_bridge
        topology = get_human_cancer_topology()
        bridge = get_human_cancer_bridge()
        config.model_name = "recon1" # Standard human model
        print("[*] Human Cancer Mode enabled (PI3K/AKT/mTOR/AMPK + Recon1)")
    elif diabetes:
        from drift.presets import get_diabetes_topology, get_diabetes_bridge
        topology = get_diabetes_topology()
        bridge = get_diabetes_bridge()
        config.model_name = "recon1"
        print("[*] Type 2 Diabetes Mode enabled (Insulin -> PI3K -> AS160 + Recon1)")

    print(f"[*] Configuration: Kd={config.drug_kd}, [Drug]={config.drug_concentration}")
    print(f"[*] Units: {config.concentration_unit}, {config.time_unit}")
    print(f"[*] Target Model: {config.model_name}")

    try:
        # Initialize Workbench
        workbench = Workbench(
            drug_kd=config.drug_kd,
            drug_concentration=config.drug_concentration,
            model_name=config.model_name,
            topology=topology,
            bridge=bridge,
            time_unit=config.time_unit,
            concentration_unit=config.concentration_unit
        )

        if run_sensitivity:
            results = workbench.run_global_sensitivity_analysis(
                n_sims=config.mc_iterations, steps=config.sim_steps, n_jobs=config.n_jobs
            )
        else:
            results = workbench.run_monte_carlo(
                n_sims=config.mc_iterations, steps=config.sim_steps, n_jobs=config.n_jobs
            )
            
        all_histories = results["histories"]
        basal_growth = results["basal_growth"]

        print(f"[+] Completed {len(all_histories)} simulations")

        # Export JSON if requested
        if json_export:
            workbench.export_results_json(results, json_export)
            print(f"[*] Results exported to JSON: {json_export}")

        # Calculate final mean growth to show in CLI
        final_growths = [h["growth"][-1] for h in all_histories]
        avg_final_growth = np.mean(final_growths)
        std_final_growth = np.std(final_growths)
        inhibition = all_histories[0]["inhibition"] * 100

        # Calculate relative growth (normalization)
        rel_growth = (avg_final_growth / basal_growth) * 100 if basal_growth > 0 else 0

        # Phenotypic stability: Coefficient of Variation
        cv = (std_final_growth / avg_final_growth) * 100 if avg_final_growth > 0 else 0

        print("\n--- Simulation Summary ---")
        print(f"Target Inhibition:      {inhibition:.2f}%")
        print(f"Basal Growth Rate:      {basal_growth:.4f} {config.time_unit}⁻¹")
        print(f"Mean Treated Growth:    {avg_final_growth:.4f} {config.time_unit}⁻¹")
        print(f"Phenotypic Uncertainty: ±{std_final_growth:.4f} ({cv:.1f}% CV)")
        print(f"Relative Vitality:      {rel_growth:.1f}%")
        print("--------------------------\n")

        print("[*] Generating interactive multi-scale dashboard...")
        fig = create_dashboard(results)

        # Pareto: Use parameter-specific filename to allow side-by-side comparison
        output_file = config.get_dashboard_filename()
        fig.write_html(output_file)

        # Also update a 'latest' symlink or copy for convenience
        latest_file = "outputs/drift_dashboard_latest.html"
        fig.write_html(latest_file)

        print("[SUCCESS] Simulation complete.")
        print(f"Results saved to: {os.path.abspath(output_file)}")
        print(f"Latest view updated: {os.path.abspath(latest_file)}")

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
    parser.add_argument(
        "--model-name", 
        type=str, 
        help="Metabolic model name (e.g., 'textbook', 'recon1', 'iJO1366')"
    )
    parser.add_argument(
        "--output", type=str, default="drift_dashboard.html", help="Output file name"
    )
    parser.add_argument(
        "--sensitivity", action="store_true", help="Run Global Sensitivity Analysis"
    )
    parser.add_argument(
        "--human-cancer", action="store_true", help="Use human cancer signaling-metabolism bridge"
    )
    parser.add_argument(
        "--diabetes", action="store_true", help="Use Type 2 Diabetes signaling-metabolism bridge"
    )
    parser.add_argument(
        "--export-json", type=str, help="Export results to JSON file"
    )
    parser.add_argument(
        "--web-config", "--web-schema", dest="web_schema", type=str, help="Export web configuration schema to JSON file"
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

    main(
        config, 
        run_sensitivity=args.sensitivity, 
        json_export=args.export_json, 
        human_cancer=args.human_cancer,
        diabetes=args.diabetes,
        web_schema=args.web_schema
    )
