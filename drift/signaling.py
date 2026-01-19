import numpy as np
from numba import njit


from .topology import Topology, get_default_topology


@njit
def default_drift_fn(state, params, feedback=None):
    """
    Default drift function for PI3K/AKT/mTOR with metabolic feedback.
    state: [PI3K, AKT, mTOR]
    params: [k_pi3k_base, k_pi3k_deg, k_akt_act, k_akt_deact, k_mtor_act, k_mtor_deact, inhibition]
    feedback: Vector of feedback signals (0-1), one for each species. 
              mTOR (index 2) is the primary target for metabolic feedback.
    """
    pi3k, akt, mtor = state
    (
        k_pi3k_base,
        k_pi3k_deg,
        k_akt_act,
        k_akt_deact,
        k_mtor_act,
        k_mtor_deact,
        inhibition,
    ) = params

    # Default feedback to 1.0 if not provided
    # Numba requires explicit types or initialization for arrays
    fb = np.ones(3)
    if feedback is not None:
        fb = feedback

    # Drug inhibits PI3K activity/availability
    effective_pi3k = pi3k * (1.0 - inhibition)

    # PI3K dynamics (basal synthesis and degradation)
    dpi3k = k_pi3k_base - k_pi3k_deg * pi3k

    # AKT dynamics (activated by PI3K)
    dakt = k_akt_act * effective_pi3k * (1.0 - akt) - k_akt_deact * akt

    # mTOR dynamics (activated by AKT, modulated by metabolic feedback)
    # Primary metabolic feedback acts on mTOR
    effective_mtor_act = k_mtor_act * fb[2]
    dmtor = effective_mtor_act * akt * (1.0 - mtor) - k_mtor_deact * mtor

    return np.array([dpi3k, dakt, dmtor])


@njit
def langevin_step(state, dt, params, noise_scale, feedback=None):
    """
    One step of SDE integration using the Milstein scheme for better stability.
    Uses state-dependent noise to ensure biological plausibility near [0, 1] boundaries.
    """
    drift = default_drift_fn(state, params, feedback=feedback)
    
    # Random term
    dw = np.random.normal(0, 1.0, size=len(state)) * np.sqrt(dt)
    
    # State-dependent noise: b(x) = noise_scale * sqrt(x * (1-x))
    # This naturally vanishes at boundaries, improving stability.
    # We use a small epsilon to avoid sqrt(0) and division by zero.
    eps = 1e-6
    clamped_state = np.zeros_like(state)
    for i in range(len(state)):
        clamped_state[i] = max(eps, min(1.0 - eps, state[i]))
        
    b = noise_scale * np.sqrt(clamped_state * (1.0 - clamped_state))
    
    # Milstein term: 0.5 * b(x) * b'(x) * (dw^2 - dt)
    # b'(x) = noise_scale * (1 - 2x) / (2 * sqrt(x * (1 - x)))
    # b(x) * b'(x) = noise_scale^2 * (1 - 2x) / 2
    milstein_corr = 0.5 * (noise_scale**2) * (1.0 - 2.0 * clamped_state) * (dw**2 - dt)

    # Update state
    new_state = state + drift * dt + b * dw + milstein_corr

    # Physical constraints: proteins stay in [0, 1] range (normalized)
    for i in range(len(new_state)):
        if new_state[i] < 0:
            new_state[i] = 0
        if new_state[i] > 1:
            new_state[i] = 1

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
        milstein_corr = 0.5 * (noise_scale**2) * (1.0 - 2.0 * clamped_state) * (dw**2 - dt)

        # Update state
        new_state = state + drift * dt + b * dw + milstein_corr

        # Physical constraints: proteins stay in [0, 1] range
        for i in range(len(new_state)):
            if new_state[i] < 0:
                new_state[i] = 0
            if new_state[i] > 1:
                new_state[i] = 1

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

        # Performance optimization: if using default topology, use jitted function
        if self.topology.name == "PI3K_AKT_mTOR":
            params = np.concatenate((self.base_params, [inhibition]))
            return langevin_step(state, self.dt, params, self.noise_scale, feedback=feedback)

        # High-performance custom topology
        if self.topology.jitted_step_fn is not None:
            # We assume the jitted_step_fn has signature (state, dt, params, noise_scale, feedback)
            # or it handles inhibition inside its params
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
