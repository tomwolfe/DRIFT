"""
Centralized audit logging for simulation decisions, SDE step-size adjustments, and constraint violations.
"""
import logging
from datetime import datetime
from typing import List, Dict, Any, Optional
from dataclasses import dataclass
from enum import Enum


class AuditSeverity(Enum):
    """Enumeration for audit log severity levels."""
    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


@dataclass
class AuditEvent:
    """Represents a single audit event."""
    timestamp: datetime
    step: int
    component: str
    message: str
    severity: AuditSeverity
    metadata: Optional[Dict[str, Any]] = None


class SimulationAuditLog:
    """
    Centralized audit logging system for capturing solver decisions, SDE step-size adjustments,
    and constraint violations for every simulation step.
    """
    
    def __init__(self):
        self.events: List[AuditEvent] = []
        self.logger = logging.getLogger(__name__)
    
    def log_event(
        self, 
        step: int, 
        component: str, 
        message: str, 
        severity: AuditSeverity = AuditSeverity.INFO,
        metadata: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        Log an audit event with timestamp, step, component, message, severity, and optional metadata.
        
        Args:
            step: Current simulation step number
            component: Component that generated the event (e.g., 'FBA Solver', 'SDE Integrator')
            message: Description of the event
            severity: Severity level of the event
            metadata: Optional additional data associated with the event
        """
        event = AuditEvent(
            timestamp=datetime.now(),
            step=step,
            component=component,
            message=message,
            severity=severity,
            metadata=metadata or {}
        )
        self.events.append(event)
        
        # Also log to standard logging for debugging
        log_method = getattr(self.logger, severity.value)
        log_method(f"[Step {step}] {component}: {message}")
    
    def log_solver_decision(
        self, 
        step: int, 
        status: str, 
        objective_value: float, 
        diagnostic: Optional[str] = None
    ) -> None:
        """Log FBA solver decision with status and objective value."""
        message = f"FBA solver status: {status}, objective: {objective_value:.4f}"
        if diagnostic:
            message += f", diagnostic: {diagnostic}"
            
        severity = AuditSeverity.WARNING if status != "optimal" else AuditSeverity.INFO
        self.log_event(
            step=step,
            component="FBA Solver",
            message=message,
            severity=severity,
            metadata={
                "status": status,
                "objective_value": objective_value,
                "diagnostic": diagnostic
            }
        )
    
    def log_constraint_violation(
        self, 
        step: int, 
        constraint_id: str, 
        expected_value: float, 
        actual_value: float
    ) -> None:
        """Log constraint violation with expected vs actual values."""
        message = f"Constraint violation: {constraint_id}, expected: {expected_value}, actual: {actual_value}"
        self.log_event(
            step=step,
            component="Constraint Validator",
            message=message,
            severity=AuditSeverity.WARNING,
            metadata={
                "constraint_id": constraint_id,
                "expected_value": expected_value,
                "actual_value": actual_value,
                "violation_amount": abs(expected_value - actual_value)
            }
        )
    
    def log_step_size_adjustment(
        self, 
        step: int, 
        old_step_size: float, 
        new_step_size: float, 
        reason: str
    ) -> None:
        """Log SDE step-size adjustment."""
        message = f"Step size adjusted: {old_step_size} -> {new_step_size}, reason: {reason}"
        self.log_event(
            step=step,
            component="SDE Integrator",
            message=message,
            severity=AuditSeverity.INFO,
            metadata={
                "old_step_size": old_step_size,
                "new_step_size": new_step_size,
                "reason": reason
            }
        )
    
    def get_events_by_severity(self, severity: AuditSeverity) -> List[AuditEvent]:
        """Get all events of a specific severity level."""
        return [event for event in self.events if event.severity == severity]
    
    def get_events_by_component(self, component: str) -> List[AuditEvent]:
        """Get all events from a specific component."""
        return [event for event in self.events if event.component == component]
    
    def get_summary_stats(self) -> Dict[str, int]:
        """Get summary statistics of audit events."""
        stats = {
            "total_events": len(self.events),
            "debug_count": 0,
            "info_count": 0,
            "warning_count": 0,
            "error_count": 0,
            "critical_count": 0
        }
        
        for event in self.events:
            stats[f"{event.severity.value}_count"] += 1
            
        return stats
    
    def clear(self) -> None:
        """Clear all audit events."""
        self.events.clear()