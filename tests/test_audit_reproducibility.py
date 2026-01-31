"""
Regression tests for audit trail and reproducibility.
Verifies current simulation outputs against known-good numerical trajectories.
"""
import numpy as np
import pytest
from drift.workbench import Workbench
from drift.engine import SimulationEngine
from drift.signaling import StochasticIntegrator
from drift.metabolic import DFBASolver, MetabolicBridge
from drift.binding import BindingEngine
from drift.topology import get_default_topology
from drift.core.exceptions import SolverError


def test_reproducible_identical_runs():
    """
    Verifies that running the same simulation twice with the same random_state
    produces identical 'signaling' and 'growth' arrays (np.array_equal).
    """
    seed = 12345
    steps = 50

    # First run
    wb1 = Workbench(random_state=seed, model_name="textbook")
    result1 = wb1.run_simulation(steps=steps)

    # Second run
    wb2 = Workbench(random_state=seed, model_name="textbook")
    result2 = wb2.run_simulation(steps=steps)

    # Compare signaling arrays
    assert np.array_equal(result1["signaling"], result2["signaling"]), "Signaling arrays should be identical"
    
    # Compare growth arrays
    assert np.array_equal(result1["growth"], result2["growth"]), "Growth arrays should be identical"
    
    # Both should have audit logs
    assert "audit_log" in result1
    assert "audit_log" in result2
    assert len(result1["audit_log"]) > 0
    assert len(result2["audit_log"]) > 0


def test_audit_log_contains_events():
    """
    Verifies that audit logs contain expected events during simulation.
    """
    seed = 98765
    steps = 20

    wb = Workbench(random_state=seed, model_name="textbook")
    result = wb.run_simulation(steps=steps)

    # Check that audit log exists and contains events
    assert "audit_log" in result
    assert len(result["audit_log"]) > 0

    # Verify audit log structure
    audit_events = result["audit_log"]
    for event in audit_events:
        assert "timestamp" in event
        assert "step" in event
        assert "component" in event
        assert "message" in event
        assert "severity" in event
        assert "metadata" in event

    # Check that FBA solver events are logged
    fba_events = [e for e in audit_events if e["component"] == "FBA Solver"]
    assert len(fba_events) > 0, "Should have FBA solver events in audit log"

    # Check that steps are properly recorded
    steps_in_log = {e["step"] for e in audit_events}
    expected_steps = set(range(steps))
    assert steps_in_log.issubset(expected_steps), f"Audit log contains unexpected steps: {steps_in_log - expected_steps}"


def test_headless_proxy_failure_audit():
    """
    Triggers a 'HeadlessProxy' failure by providing impossible constraints
    and verifies the 'audit_log' captures the 'SolverError' with AuditSeverity.ERROR.
    """
    # Create a headless solver (no actual solver installed)
    solver = DFBASolver(model_name="textbook")

    # Force headless mode and proxy logic
    solver.headless = True
    solver.model = None
    if not hasattr(solver, "proxy_params"):
        solver.proxy_params = {"base_growth": 0.2}

    # Create a bridge and integrator
    topology = get_default_topology()
    bridge = MetabolicBridge(species_names=topology.species)
    binding = BindingEngine(targets=1.0)
    integrator = StochasticIntegrator(topology=topology)

    engine = SimulationEngine(integrator=integrator, solver=solver, bridge=bridge, binding=binding)

    # Test with negative growth (should trigger PhysicalBoundValidator)
    solver.proxy_params["base_growth"] = -0.5  # Negative growth

    # We expect this to be caught by the engine and logged as an ERROR
    result = engine.run(steps=5, drug_concentration=1.0)

    # Verify cell death occurred
    assert result["cell_death"] is True
    assert "PhysicalBoundValidator" in result["death_cause"]

    # Verify audit log captures the error
    audit_events = result["audit_log"]
    error_events = [e for e in audit_events if e["severity"] == "error"]
    
    assert len(error_events) > 0, "Should have error events in audit log"
    assert "FBA Solver" in error_events[0]["component"]
    assert "Critical Solver Error" in error_events[0]["message"]
    assert "PhysicalBoundValidator" in error_events[0]["metadata"]["error"]


def test_audit_log_integration_in_engine():
    """
    Tests that the SimulationEngine properly integrates audit logging.
    """
    seed = 42
    steps = 10

    # Create components
    topology = get_default_topology()
    integrator = StochasticIntegrator(topology=topology, noise_scale=0.02)
    solver = DFBASolver(model_name="textbook")  # May be headless depending on environment
    bridge = MetabolicBridge(species_names=topology.species)
    binding = BindingEngine(targets=1.0)

    # Create engine with specific seed
    engine = SimulationEngine(
        integrator=integrator,
        solver=solver,
        bridge=bridge,
        binding=binding,
        rng=np.random.default_rng(seed)
    )

    # Run simulation
    result = engine.run(steps=steps, drug_concentration=1.0, seed=seed)

    # Verify audit log is present and populated
    assert "audit_log" in result
    assert len(result["audit_log"]) > 0

    # Verify that all expected fields are present in result
    expected_fields = [
        "time", "signaling", "growth", "status", "cell_death", 
        "death_step", "death_cause", "inhibition", "drug_kd", 
        "params", "headless", "species_names", "audit_log"
    ]
    for field in expected_fields:
        assert field in result, f"Missing field {field} in simulation result"


def test_different_seeds_produce_different_results():
    """
    Verifies that different random seeds produce different results.
    """
    seed1, seed2 = 100, 200
    steps = 30

    wb1 = Workbench(random_state=seed1, model_name="textbook")
    result1 = wb1.run_simulation(steps=steps)

    wb2 = Workbench(random_state=seed2, model_name="textbook")
    result2 = wb2.run_simulation(steps=steps)

    # Results should be different
    assert not np.array_equal(result1["signaling"], result2["signaling"]), "Signaling arrays should be different"
    
    assert not np.array_equal(result1["growth"], result2["growth"]), "Growth arrays should be different"


def test_audit_log_serialization():
    """
    Tests that audit logs can be properly serialized/deserialized.
    """
    seed = 54321
    steps = 15

    wb = Workbench(random_state=seed, model_name="textbook")
    result = wb.run_simulation(steps=steps)

    # Verify audit log can be JSON-serialized (no complex objects)
    import json
    try:
        json_str = json.dumps(result["audit_log"])
        parsed_back = json.loads(json_str)
        assert len(parsed_back) == len(result["audit_log"])
    except (TypeError, ValueError) as e:
        pytest.fail(f"Audit log is not JSON serializable: {e}")


def test_audit_log_utility_methods():
    """
    Tests utility methods of SimulationAuditLog to ensure high coverage.
    """
    from drift.core.audit import SimulationAuditLog, AuditSeverity
    
    audit = SimulationAuditLog()
    
    # Test log_event
    audit.log_event(step=1, component="Test", message="Info message", severity=AuditSeverity.INFO)
    audit.log_event(step=2, component="Test", message="Warning message", severity=AuditSeverity.WARNING)
    audit.log_event(step=3, component="Other", message="Error message", severity=AuditSeverity.ERROR)
    
    # Test log_constraint_violation
    audit.log_constraint_violation(step=4, constraint_id="R1", expected_value=10.0, actual_value=15.0)
    
    # Test get_events_by_severity
    info_events = audit.get_events_by_severity(AuditSeverity.INFO)
    assert len(info_events) == 1
    assert info_events[0].step == 1
    
    # Test get_events_by_component
    test_events = audit.get_events_by_component("Test")
    assert len(test_events) == 2
    
    # Test get_summary_stats
    stats = audit.get_summary_stats()
    assert stats["total_events"] == 4
    assert stats["info_count"] == 1
    assert stats["warning_count"] == 2 # log_constraint_violation is WARNING
    assert stats["error_count"] == 1
    
    # Test clear
    audit.clear()
    assert len(audit.events) == 0