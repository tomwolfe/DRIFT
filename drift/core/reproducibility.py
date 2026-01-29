import numpy as np
from typing import Optional

def get_rng(seed: Optional[int] = None) -> np.random.Generator:
    """
    Centralized factory for creating NumPy random generators.
    Ensures thread-safe reproducibility across simulation engines.
    """
    return np.random.default_rng(seed)
