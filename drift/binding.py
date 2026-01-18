import numpy as np

class BindingEngine:
    """Calculates protein inhibition percentage based on Kd and [Drug]."""
    def __init__(self, kd: float, hill_coefficient: float = 1.0):
        self.kd = kd
        self.hill_coefficient = hill_coefficient

    def calculate_occupancy(self, drug_concentration: float) -> float:
        """Returns the fraction of target bound by drug [0, 1]."""
        if drug_concentration <= 0:
            return 0.0
        
        # Hill Equation: theta = [D]^n / ([D]^n + Kd^n)
        return (drug_concentration ** self.hill_coefficient) / \
               (drug_concentration ** self.hill_coefficient + self.kd ** self.hill_coefficient)

    def calculate_inhibition(self, drug_concentration: float) -> float:
        """Returns inhibition percentage [0.0, 1.0]."""
        return self.calculate_occupancy(drug_concentration)
