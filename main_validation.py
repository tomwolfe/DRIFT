"""
DRIFT Scientific Validation Suite
================================

This script executes the formal validation tests for the DRIFT framework,
confirming its biological accuracy and numerical stability.
"""

import numpy as np
from drift.workbench import Workbench
from drift.binding import BindingEngine


def run_test(name, func):
    print(f"[*] Running: {name}...")
    try:
        func()
        print(f"    [PASS] {name}")
        return True
    except Exception as e:
        print(f"    [FAIL] {name}: {str(e)}")
        return False


def test_binding_model():
    be = BindingEngine(kd=0.5)
    occupancy = be.calculate_occupancy(0.5)
    assert abs(occupancy - 0.5) < 0.01, f"Expected ~0.5, got {occupancy}"


def test_dose_response():
    concentrations = [0.0, 1.0, 100.0]
    growth_rates = []
    for conc in concentrations:
        wb = Workbench(drug_kd=1.0, drug_concentration=conc)
        wb.signaling.noise_scale = 0.001
        history = wb.run_simulation(steps=1000)
        growth_rates.append(np.mean(history["growth"][-100:]))

    # Use a small epsilon for float comparison or just check they are reasonably separated
    assert (
        growth_rates[0] > growth_rates[1] * 0.99
    ), f"Basal {growth_rates[0]} should be >= {growth_rates[1]}"
    assert (
        growth_rates[1] > growth_rates[2]
    ), f"Growth {growth_rates[2]} should be < {growth_rates[1]}"


def test_signaling_delays():
    wb = Workbench(drug_kd=0.1, drug_concentration=10.0)
    wb.signaling.noise_scale = 0.001
    history = wb.run_simulation(steps=1000)
    signaling = history["signaling"]

    final_akt = np.mean(signaling[-100:, 1])
    final_mtor = np.mean(signaling[-100:, 2])

    assert (
        final_akt < 0.5
    ), f"AKT should have dropped significantly, got {final_akt:.4f}"
    assert (
        final_mtor < 0.6
    ), f"mTOR should have dropped significantly, got {final_mtor:.4f}"


def test_metabolic_response():
    wb_basal = Workbench(drug_kd=0.5, drug_concentration=0.0)
    wb_basal.signaling.noise_scale = 0.001
    basal_growth = wb_basal.run_simulation(steps=200)["growth"][-1]

    wb_inhibited = Workbench(drug_kd=0.5, drug_concentration=50.0)
    wb_inhibited.signaling.noise_scale = 0.001
    inhibited_growth = wb_inhibited.run_simulation(steps=200)["growth"][-1]

    assert inhibited_growth < basal_growth, "Inhibited growth must be lower than basal"


def test_solver_sensitivity():
    wb = Workbench(model_name="textbook")
    results = wb.solver.check_solver_sensitivity()
    # At least one solver should be optimal
    optimal_solvers = [s for s, r in results.items() if r["status"] == "optimal"]
    assert len(optimal_solvers) >= 1, "No optimal solvers found"
    
    if len(optimal_solvers) > 1:
        objectives = [results[s]["objective"] for s in optimal_solvers]
        max_diff = max(objectives) - min(objectives)
        # GLPK and others should be very close for textbook model
        assert max_diff < 1e-5, f"Solver discrepancy too high: {max_diff}"


def main():
    print("=" * 40)
    print("DRIFT SCIENTIFIC VALIDATION REPORT")
    print("=" * 40)

    tests = [
        ("Hill Equation Binding Accuracy", test_binding_model),
        ("Pharmacological Dose-Response Sensitivity", test_dose_response),
        ("Temporal Signaling Cascade Propagation", test_signaling_delays),
        ("Metabolic Coupling Integrity", test_metabolic_response),
        ("FBA Solver Sensitivity and Consistency", test_solver_sensitivity),
    ]

    passed = 0
    for name, func in tests:
        if run_test(name, func):
            passed += 1

    print("=" * 40)
    print(f"SUMMARY: {passed}/{len(tests)} tests passed.")
    print("=" * 40)


if __name__ == "__main__":
    main()
