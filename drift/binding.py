class BindingEngine:
    """Calculates protein inhibition percentage based on Kd and [Drug]."""

    kd: float
    hill_coefficient: float

    def __init__(self, kd: float, hill_coefficient: float = 1.0):
        """
        Initialize the BindingEngine.

        Args:
            kd (float): Dissociation constant of the drug
            hill_coefficient (float): Hill coefficient for cooperativity (default 1.0)

        Raises:
            ValueError: If kd is not positive
        """
        if kd <= 0:
            raise ValueError(f"kd must be positive, got {kd}")
        if hill_coefficient <= 0:
            raise ValueError(
                f"hill_coefficient must be positive, got {hill_coefficient}"
            )

        self.kd = kd
        self.hill_coefficient = hill_coefficient

    def calculate_occupancy(self, drug_concentration: float) -> float:
        """
        Returns the fraction of target bound by drug [0, 1].

        Args:
            drug_concentration (float): Concentration of the drug

        Returns:
            float: Fraction of target bound by drug [0, 1]
        """
        if drug_concentration < 0:
            raise ValueError(
                f"drug_concentration must be non-negative, got {drug_concentration}"
            )

        if drug_concentration == 0:
            return 0.0

        # Hill Equation: theta = [D]^n / ([D]^n + Kd^n)
        numerator = drug_concentration**self.hill_coefficient
        denominator = numerator + self.kd**self.hill_coefficient
        return float(numerator / denominator)

    def calculate_inhibition(self, drug_concentration: float) -> float:
        """
        Returns inhibition percentage [0.0, 1.0].

        Args:
            drug_concentration (float): Concentration of the drug

        Returns:
            float: Inhibition percentage [0.0, 1.0]
        """
        return float(self.calculate_occupancy(drug_concentration))
