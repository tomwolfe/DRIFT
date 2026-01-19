import numpy as np
from numba import njit


from .topology import Topology, get_default_topology


@njit
def default_drift_fn(state, params):
    """
    Default drift function for PI3K/AKT/mTOR.
    state: [PI3K, AKT, mTOR]
    params: [k_pi3k_base, k_pi3k_deg, k_akt_act, k_akt_deact, k_mtor_act, k_mtor_deact, inhibition]
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

    # Drug inhibits PI3K activity/availability
    effective_pi3k = pi3k * (1.0 - inhibition)

    # PI3K dynamics (basal synthesis and degradation)
    dpi3k = k_pi3k_base - k_pi3k_deg * pi3k

    # AKT dynamics (activated by PI3K)
    dakt = k_akt_act * effective_pi3k * (1.0 - akt) - k_akt_deact * akt

    # mTOR dynamics (activated by AKT)
    dmtor = k_mtor_act * akt * (1.0 - mtor) - k_mtor_deact * mtor

    return np.array([dpi3k, dakt, dmtor])


@njit
def langevin_step(state, dt, params, noise_scale):
    """
    One step of Langevin dynamics using the default drift function.
    """
    drift = default_drift_fn(state, params)

    # Update state
    new_state = state + drift * dt + np.random.normal(0, noise_scale, size=len(state)) * np.sqrt(dt)

    # Physical constraints: proteins stay in [0, 1] range (normalized)
    for i in range(len(new_state)):
        if new_state[i] < 0:
            new_state[i] = 0
        if new_state[i] > 1:
            new_state[i] = 1

    return new_state


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

    def step(self, state, inhibition):
        """
        Perform one integration step.
        """
        if not isinstance(state, np.ndarray) or state.shape != (len(self.topology.species),):
            raise ValueError(
                f"state must be a numpy array of shape ({len(self.topology.species)},), got {type(state)} with shape {getattr(state, 'shape', 'N/A')}"
            )
        
        # Performance optimization: if using default topology, use jitted function
        if self.topology.name == "PI3K_AKT_mTOR":
            params = np.concatenate((self.base_params, [inhibition]))
            return langevin_step(state, self.dt, params, self.noise_scale)
        
        # Custom topology with provided drift_fn
        if self.topology.drift_fn is not None:
            return self._custom_step(state, inhibition)
            
        # Fallback to generic decay
        return self._generic_step(state, inhibition)

    def _custom_step(self, state, inhibition):
        """Step using a custom drift function provided in the topology."""
        drift = self.topology.drift_fn(state, self.topology.parameters, inhibition)
        diffusion = np.random.normal(0, self.noise_scale, size=len(state)) * np.sqrt(self.dt)
        new_state = state + drift * self.dt + diffusion
        return np.clip(new_state, 0, 1)

    def _generic_step(self, state, inhibition):
        """Generic Euler-Maruyama step for unknown topologies."""
        # Fallback: very simple decay-to-zero drift
        drift = -0.1 * state
        diffusion = np.random.normal(0, self.noise_scale, size=len(state)) * np.sqrt(self.dt)
        new_state = state + drift * self.dt + diffusion
        return np.clip(new_state, 0, 1)
