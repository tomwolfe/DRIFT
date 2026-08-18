#!/usr/bin/env python
"""
DRIFT Performance Benchmark Suite

Runs systematic performance benchmarks across different configurations
and reports timing, memory usage, and solver performance.
"""

import time
import os
import sys
import csv
import gc
from pathlib import Path
from typing import Dict, Any, List, Optional

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from drift.workbench import Workbench
from drift.binding import BindingEngine
from drift.topology import get_default_topology


def get_memory_usage_mb() -> float:
    """Get current process memory usage in MB."""
    try:
        import resource
        # macOS/Linux: ru_maxrss is in bytes on macOS, KB on Linux
        usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform == "darwin":
            return usage / (1024 * 1024)  # bytes to MB
        else:
            return usage / 1024  # KB to MB
    except ImportError:
        return 0.0


def benchmark_single_run(
    steps: int,
    model_name: str = "textbook",
    drug_kd: float = 1.0,
    drug_concentration: float = 2.0,
    sub_steps: int = 5,
    seed: int = 42,
) -> Dict[str, Any]:
    """Benchmark a single simulation run."""
    gc.collect()
    mem_before = get_memory_usage_mb()

    wb = Workbench(
        drug_kd=drug_kd,
        drug_concentration=drug_concentration,
        model_name=model_name,
        sub_steps=sub_steps,
    )
    wb.signaling.noise_scale = 0.001

    start = time.perf_counter()
    history = wb.run_simulation(steps=steps, seed=seed)
    elapsed = time.perf_counter() - start

    mem_after = get_memory_usage_mb()

    return {
        "model": model_name,
        "steps": steps,
        "sub_steps": sub_steps,
        "n_sims": 1,
        "wall_time_s": elapsed,
        "time_per_step_ms": (elapsed / steps) * 1000,
        "memory_mb": mem_after - mem_before,
        "final_growth": float(history["growth"][-1]),
        "cell_death": history.get("cell_death", False),
    }


def benchmark_monte_carlo(
    n_sims: int,
    steps: int,
    model_name: str = "textbook",
    n_jobs: int = 1,
    seed: int = 42,
) -> Dict[str, Any]:
    """Benchmark a Monte Carlo ensemble run."""
    gc.collect()
    mem_before = get_memory_usage_mb()

    wb = Workbench(
        drug_kd=1.0,
        drug_concentration=2.0,
        model_name=model_name,
        sub_steps=5,
        random_state=seed,
    )
    wb.signaling.noise_scale = 0.001

    start = time.perf_counter()
    results = wb.run_monte_carlo(n_sims=n_sims, steps=steps, n_jobs=n_jobs)
    elapsed = time.perf_counter() - start

    mem_after = get_memory_usage_mb()

    if isinstance(results, dict) and "histories" in results:
        n_deaths = sum(1 for h in results["histories"] if h.get("cell_death", False))
    else:
        n_deaths = 0

    return {
        "model": model_name,
        "steps": steps,
        "n_sims": n_sims,
        "n_jobs": n_jobs,
        "wall_time_s": elapsed,
        "time_per_sim_s": elapsed / n_sims,
        "memory_mb": mem_after - mem_before,
        "n_deaths": n_deaths,
    }


def run_benchmarks(output_path: Optional[str] = None) -> List[Dict[str, Any]]:
    """Run all benchmarks and save results."""
    if output_path is None:
        output_path = str(Path(__file__).parent.parent / "outputs" / "benchmark_results.csv")

    results = []

    print("=" * 60)
    print("DRIFT Performance Benchmark Suite")
    print("=" * 60)

    # 1. Single-run benchmarks with different step counts
    print("\n--- Single-Run Benchmarks (textbook model) ---")
    for steps in [10, 50, 100]:
        print(f"  Steps={steps}...", end=" ", flush=True)
        result = benchmark_single_run(steps=steps, model_name="textbook")
        results.append(result)
        print(f"{result['wall_time_s']:.2f}s ({result['time_per_step_ms']:.1f} ms/step)")

    # 2. Monte Carlo benchmarks with different ensemble sizes
    print("\n--- Monte Carlo Benchmarks (textbook model) ---")
    for n_sims in [5, 10, 30]:
        print(f"  N={n_sims}, steps=50...", end=" ", flush=True)
        result = benchmark_monte_carlo(n_sims=n_sims, steps=50, n_jobs=1)
        results.append(result)
        print(f"{result['wall_time_s']:.2f}s ({result['time_per_sim_s']:.2f} s/sim)")

    # 3. Parallel benchmark
    print("\n--- Parallel Benchmark (textbook, N=10) ---")
    for n_jobs in [1, 2, 4]:
        n_cpus = os.cpu_count() or 1
        if n_jobs > n_cpus:
            print(f"  Skipping n_jobs={n_jobs} (only {n_cpus} CPUs available)")
            continue
        print(f"  n_jobs={n_jobs}...", end=" ", flush=True)
        result = benchmark_monte_carlo(n_sims=10, steps=50, n_jobs=n_jobs)
        results.append(result)
        print(f"{result['wall_time_s']:.2f}s")

    # Save results
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    if results:
        # Collect all unique keys across all results
        all_keys = set()
        for r in results:
            all_keys.update(r.keys())
        all_keys = sorted(all_keys)

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_keys)
            writer.writeheader()
            for r in results:
                writer.writerow({k: r.get(k, "") for k in all_keys})
        print(f"\n[+] Results saved to {output_path}")

    # Summary
    print("\n--- Summary ---")
    print(f"Total benchmarks: {len(results)}")
    if results:
        single_runs = [r for r in results if r.get("n_sims", 0) == 1]
        if single_runs:
            print(f"Single-run time range: {min(r['wall_time_s'] for r in single_runs):.2f}s - {max(r['wall_time_s'] for r in single_runs):.2f}s")
        mc_runs = [r for r in results if r.get("n_sims", 0) > 1]
        if mc_runs:
            print(f"MC ensemble time range: {min(r['wall_time_s'] for r in mc_runs):.2f}s - {max(r['wall_time_s'] for r in mc_runs):.2f}s")

    return results


if __name__ == "__main__":
    run_benchmarks()
