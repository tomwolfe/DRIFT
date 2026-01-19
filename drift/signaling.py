import numpy as np
from numba import njit


from .topology import Topology, get_default_topology


@njit
def langevin_step(state, dt, params, noise_scale):
    """
    One step of Langevin dynamics for PI3K/AKT/mTOR.
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

    # Update state
    drift = np.array([dpi3k, dakt, dmtor]) * dt
    diffusion = np.random.normal(0, noise_scale, size=3) * np.sqrt(dt)

    new_state = state + drift + diffusion

    # Physical constraints: proteins stay in [0, 1] range (normalized)
    for i in range(len(state)):
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

        Raises:
            ValueError: If dt or noise_scale are invalid
        """
        if dt <= 0:
            raise ValueError(f"dt must be positive, got {dt}")
        if noise_scale < 0:
            raise ValueError(f"noise_scale must be non-negative, got {noise_scale}")

        self.dt = dt
        self.noise_scale = noise_scale
        self.topology = topology or get_default_topology()
        
        # Default kinetic parameters for PI3K/AKT/mTOR
        self.base_params = np.array(list(self.topology.parameters.values()))

    def step(self, state, inhibition):
        """
        Perform one integration step.

        Args:
            state (np.ndarray): Current state
            inhibition (float): Inhibition factor from binding

        Returns:
            np.ndarray: Updated state after one step
        """
        if not isinstance(state, np.ndarray) or state.shape != (len(self.topology.species),):
            raise ValueError(
                f"state must be a numpy array of shape ({len(self.topology.species)},), got {type(state)} with shape {getattr(state, 'shape', 'N/A')}"
            )
        if not (0 <= inhibition <= 1):
            raise ValueError(f"inhibition must be between 0 and 1, got {inhibition}")

        # If it's the default topology, use the fast jitted step
        if self.topology.name == "PI3K_AKT_mTOR":
            params = np.concatenate((self.base_params, [inhibition]))
            return langevin_step(state, self.dt, params, self.noise_scale)
        
        # Otherwise use a generic step (can be extended for custom ODEs)
        return self._generic_step(state, inhibition)

    def _generic_step(self, state, inhibition):
        """Generic Euler-Maruyama step for custom topologies."""
        # Simple placeholder for custom drift if provided, otherwise assume identity
        if hasattr(self.topology, 'drift_fn') and self.topology.drift_fn:
            drift = self.topology.drift_fn(state, self.topology.parameters, inhibition)
        else:
            # Fallback: very simple decay-to-zero drift if no function provided
            drift = -0.1 * state
            
        diffusion = np.random.normal(0, self.noise_scale, size=len(state)) * np.sqrt(self.dt)
        new_state = state + drift * self.dt + diffusion
        
        return np.clip(new_state, 0, 1)
