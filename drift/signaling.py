import numpy as np
from numba import njit
import logging


from .topology import Topology, get_default_topology

logger = logging.getLogger(__name__)


@njit
def seed_numba(seed):
    """Sets the seed for Numba's random number generator."""
    np.random.seed(seed)


@njit
def langevin_step_generic(state, dt, params, noise_scale, drift_fn, feedback=None):
    """
    One step of SDE integration using the Milstein scheme for better stability.
    Uses state-dependent noise to ensure biological plausibility near [0, 1] boundaries.
    Accepts an arbitrary jitted drift_fn.
    """
    drift = drift_fn(state, params, feedback=feedback)

    # Random term
    dw = np.random.normal(0, 1.0, size=len(state)) * np.sqrt(dt)

    # State-dependent noise: b(x) = noise_scale * sqrt(x * (1-x))
    eps = 1e-6
    clamped_state = np.zeros_like(state)
    for i in range(len(state)):
        clamped_state[i] = max(eps, min(1.0 - eps, state[i]))

    b = noise_scale * np.sqrt(clamped_state * (1.0 - clamped_state))
    milstein_corr = 0.25 * (noise_scale**2) * (1.0 - 2.0 * clamped_state) * (dw**2 - dt)

    # Update state
    new_state = state + drift * dt + b * dw + milstein_corr

    # Reflection principle
    for i in range(len(new_state)):
        if new_state[i] < 0:
            new_state[i] = -new_state[i]
        if new_state[i] > 1:
            new_state[i] = 2.0 - new_state[i]

        if new_state[i] < 0: new_state[i] = eps
        if new_state[i] > 1: new_state[i] = 1.0 - eps

    return new_state


@njit
def langevin_step(state, dt, params, noise_scale, feedback=None):
    """
    One step of SDE integration using the Milstein scheme for better stability.
    Uses state-dependent noise to ensure biological plausibility near [0, 1] boundaries.
    This is a simplified version that uses the default PI3K_AKT_mTOR drift function.
    """
    # Default PI3K_AKT_mTOR drift function
    def default_drift_fn(state, params, feedback=None):
        # Extract parameters for PI3K-AKT-mTOR pathway
        # params[0:3] are degradation rates, params[3:6] are production rates, params[6] is feedback
        k_deg = params[0:3]  # degradation rates for PI3K, AKT, mTOR
        k_prod = params[3:6]  # production rates
        feedback_param = params[6] if len(params) > 6 else 1.0

        if feedback is not None:
            k_prod = k_prod * feedback  # Apply feedback modulation to production rates

        # Simple linearized model for PI3K-AKT-mTOR pathway
        # dx/dt = k_prod * (1-x) - k_deg * x
        drift = k_prod * (1.0 - state) - k_deg * state
        return drift

    drift = default_drift_fn(state, params, feedback=feedback)

    # Random term
    dw = np.random.normal(0, 1.0, size=len(state)) * np.sqrt(dt)

    # State-dependent noise: b(x) = noise_scale * sqrt(x * (1-x))
    eps = 1e-6
    clamped_state = np.zeros_like(state)
    for i in range(len(state)):
        clamped_state[i] = max(eps, min(1.0 - eps, state[i]))

    b = noise_scale * np.sqrt(clamped_state * (1.0 - clamped_state))
    milstein_corr = 0.25 * (noise_scale**2) * (1.0 - 2.0 * clamped_state) * (dw**2 - dt)

    # Update state
    new_state = state + drift * dt + b * dw + milstein_corr

    # Reflection principle
    for i in range(len(new_state)):
        if new_state[i] < 0:
            new_state[i] = -new_state[i]
        if new_state[i] > 1:
            new_state[i] = 2.0 - new_state[i]

        if new_state[i] < 0: new_state[i] = eps
        if new_state[i] > 1: new_state[i] = 1.0 - eps

    return new_state


def create_langevin_integrator(drift_fn):
    """
    Higher-order function that creates a jitted Milstein integrator step 
    from a jitted drift function.
    """
    @njit
    def custom_milstein_step(state, dt, params, noise_scale, feedback=None):
        drift = drift_fn(state, params, feedback=feedback)
        
        dw = np.random.normal(0, 1.0, size=len(state)) * np.sqrt(dt)
        
        eps = 1e-6
        clamped_state = np.zeros_like(state)
        for i in range(len(state)):
            clamped_state[i] = max(eps, min(1.0 - eps, state[i]))
            
        b = noise_scale * np.sqrt(clamped_state * (1.0 - clamped_state))
        milstein_corr = 0.25 * (noise_scale**2) * (1.0 - 2.0 * clamped_state) * (dw**2 - dt)

        # Update state
        new_state = state + drift * dt + b * dw + milstein_corr

        # Reflection principle
        for i in range(len(new_state)):
            if new_state[i] < 0:
                new_state[i] = -new_state[i]
            if new_state[i] > 1:
                new_state[i] = 2.0 - new_state[i]
            
            if new_state[i] < 0: new_state[i] = eps
            if new_state[i] > 1: new_state[i] = 1.0 - eps

        return new_state
        
    return custom_milstein_step


class StochasticIntegrator:
    """Integrator for stochastic differential equations in signaling pathways."""

    def __init__(self, dt=0.1, noise_scale=0.02, topology=None):
        """
        Initialize the StochasticIntegrator.

        Args:
            dt (float): Time step for integration
            noise_scale (float): Scale of the noise term
            topology (Topology, optional): Signaling network topology
        """
        if dt <= 0:
            raise ValueError(f"dt must be positive, got {dt}")
        if noise_scale < 0:
            raise ValueError(f"noise_scale must be non-negative, got {noise_scale}")

        self.dt = dt
        self.noise_scale = noise_scale
        self.topology = topology or get_default_topology()
        
        # Default kinetic parameters
        self.base_params = np.array(list(self.topology.parameters.values()))
        
        # Performance optimization: Auto-jit custom topologies if a jitted drift is provided
        # This addresses the 'performance cliff' for custom topologies identified in audit.
        if self.topology.jitted_step_fn is None and self.topology.drift_fn is not None:
            # Check if it's already jitted (Numba CPUDispatcher) or marked by decorator
            is_jitted = type(self.topology.drift_fn).__name__ == "CPUDispatcher"
            is_marked = hasattr(self.topology.drift_fn, "_drift_model_name")
            
            if is_jitted or is_marked:
                logger.info(f"Auto-jitting Milstein integrator for custom topology: {self.topology.name}")
                try:
                    self.topology.jitted_step_fn = create_langevin_integrator(self.topology.drift_fn)
                except Exception as e:
                    logger.warning(f"Failed to auto-jit integrator for {self.topology.name}: {e}. Falling back to EM.")

        # Stability check
        self._check_stability()

    def set_seed(self, seed: int):
        """Sets the seed for the integrator (including Numba jitted functions)."""
        np.random.seed(seed)
        seed_numba(seed)
        logger.debug(f"Integrator seed set to {seed}")

    def _check_stability(self):
        """
        Heuristic stability check for the Milstein scheme.
        Warns if dt is too large relative to noise_scale or absolute time-step.
        """
        import warnings
        
        # 1. Absolute dt check
        if self.dt > 0.2:
            warnings.warn(
                f"Large time-step (dt={self.dt}) detected. This may lead to "
                "integration inaccuracies in highly non-linear signaling systems.",
                RuntimeWarning
            )
            
        # 2. Noise-driven stability check
        # For Milstein on [0, 1] bounded systems, noise_scale * sqrt(dt) 
        # should ideally be small enough to avoid frequent boundary hitting.
        noise_impact = self.noise_scale * np.sqrt(self.dt)
        if noise_impact > 0.1:
            warnings.warn(
                f"High noise impact detected (noise_scale * sqrt(dt) = {noise_impact:.3f}). "
                "Milstein scheme stability may be compromised. Consider reducing dt or noise_scale.",
                RuntimeWarning
            )

    def step(self, state, inhibition, feedback=None):
        """
        Perform one integration step with optional metabolic feedback.
        """
        if not isinstance(state, np.ndarray) or state.shape != (len(self.topology.species),):
            raise ValueError(
                f"state must be a numpy array of shape ({len(self.topology.species)},), got {type(state)} with shape {getattr(state, 'shape', 'N/A')}"
            )

        # Validate inhibition and feedback parameters
        if not (0 <= inhibition <= 1):
            raise ValueError(
                f"inhibition must be between 0 and 1, got {inhibition}"
            )
        
        if feedback is None:
            feedback = np.ones(len(self.topology.species))
        elif isinstance(feedback, (float, int)):
            # Convert scalar feedback to vector for backward compatibility
            val = float(feedback)
            feedback = np.ones(len(self.topology.species)) * val

        # Performance optimization: if using jitted topology, use it
        if self.topology.jitted_step_fn is not None:
            # We assume the jitted_step_fn has signature (state, dt, params, noise_scale, feedback)
            params = np.concatenate((self.base_params, [inhibition]))
            return self.topology.jitted_step_fn(state, self.dt, params, self.noise_scale, feedback)

        # Custom topology with provided drift_fn (non-jitted fallback)
        if self.topology.drift_fn is not None:
            return self._custom_step(state, inhibition, feedback=feedback)

        # Fallback to generic decay
        return self._generic_step(state, inhibition, feedback=feedback)

    def _custom_step(self, state, inhibition, feedback=None):
        """Step using a custom drift function provided in the topology."""
        if self.topology.drift_fn is None:
            # Fallback if drift_fn is None
            return self._generic_step(state, inhibition, feedback=feedback)

        if feedback is None:
            feedback = np.ones(len(state))

        # Custom drift functions should now ideally accept feedback as well
        # We try to pass it if possible, otherwise we ignore it for backward compatibility
        try:
            drift = self.topology.drift_fn(state, self.topology.parameters, inhibition, feedback=feedback)
        except TypeError:
            drift = self.topology.drift_fn(state, self.topology.parameters, inhibition)

        diffusion = np.random.normal(0, self.noise_scale, size=len(state)) * np.sqrt(self.dt)
        new_state = state + drift * self.dt + diffusion
        return np.clip(new_state, 0, 1)

    def _generic_step(self, state, inhibition, feedback=None):
        """Generic Euler-Maruyama step for unknown topologies."""
        if feedback is None:
            feedback = np.ones(len(state))
            
        # Fallback: simple decay modulated by mean feedback for generic case
        mean_feedback = np.mean(feedback)
        drift = -0.1 * state * (2.0 - mean_feedback) 
        
        # Apply inhibition if topology specifies an inhibited species
        inhibited_species = getattr(self.topology, "inhibited_species", None)
        if inhibited_species and inhibited_species in self.topology.species:
            idx = self.topology.species.index(inhibited_species)
            # Reduce the "activity" (clamped level) of the target species
            state_with_inhibition = state.copy()
            state_with_inhibition[idx] *= (1.0 - inhibition)
            # Recalculate drift with inhibition if we had a more complex generic model
            # For now, just apply it to the state used in the next step
            state = state_with_inhibition

        diffusion = np.random.normal(0, self.noise_scale, size=len(state)) * np.sqrt(self.dt)
        new_state = state + drift * self.dt + diffusion
        return np.clip(new_state, 0, 1)
