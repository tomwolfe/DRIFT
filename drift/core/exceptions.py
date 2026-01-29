"""
Custom exception classes for the DRIFT framework.
"""

class DriftError(Exception):
    """Base class for exceptions in DRIFT."""
    pass

class DesyncError(DriftError):
    """Exception raised when signaling and metabolic layers are desynchronized."""
    pass

class SolverError(DriftError):
    """Exception raised when the DFBASolver encounters an unrecoverable error."""
    pass

class ConfigurationError(DriftError):
    """Exception raised for errors in the configuration."""
    pass
