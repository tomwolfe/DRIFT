from typing import Union, Dict, List, Optional

class BindingEngine:
    """Calculates protein inhibition percentage based on Kd and [Drug]."""

    def __init__(self, targets: Optional[Union[float, Dict[str, float]]] = None, hill_coefficient: float = 1.0, kd: Optional[float] = None):
        """
        Initialize the BindingEngine.

        Args:
            targets (float or dict): If float, the Kd for the default inhibited species.
                                     If dict, mapping of species names to their respective Kd values.
            hill_coefficient (float): Hill coefficient for cooperativity (default 1.0)
            kd (float, optional): For backward compatibility, same as targets if targets is float.

        Raises:
            ValueError: If kd is not positive
        """
        if hill_coefficient <= 0:
            raise ValueError(
                f"hill_coefficient must be positive, got {hill_coefficient}"
            )
        
        eff_targets = targets if targets is not None else kd
        if eff_targets is None:
            raise ValueError("Must provide either targets or kd")

        if isinstance(eff_targets, (int, float)):
            if eff_targets <= 0:
                raise ValueError(f"kd must be positive, got {eff_targets}")
            self.targets = {"default": float(eff_targets)}
            self.kd = float(eff_targets) # Backward compatibility
        else:
            for name, val in eff_targets.items():
                if val <= 0:
                    raise ValueError(f"kd for {name} must be positive, got {val}")
            self.targets = {name: float(val) for name, val in eff_targets.items()}
            # For backward compatibility, use the first Kd if available
            self.kd = next(iter(self.targets.values())) if self.targets else 1.0

        self.hill_coefficient = hill_coefficient

    def calculate_occupancy(self, drug_concentration: float, target_name: str = "default") -> float:
        """
        Returns the fraction of target bound by drug [0, 1].

        Args:
            drug_concentration (float): Concentration of the drug
            target_name (str): Name of the target species to calculate occupancy for.

        Returns:
            float: Fraction of target bound by drug [0, 1]
        """
        if drug_concentration < 0:
            raise ValueError(
                f"drug_concentration must be non-negative, got {drug_concentration}"
            )

        if drug_concentration == 0:
            return 0.0

        kd = self.targets.get(target_name, self.kd)

        # Hill Equation: theta = [D]^n / ([D]^n + Kd^n)
        numerator = drug_concentration**self.hill_coefficient
        denominator = numerator + kd**self.hill_coefficient
        return float(numerator / denominator)

    def calculate_inhibition(self, drug_concentration: float) -> Union[float, Dict[str, float]]:
        """
        Returns inhibition percentage [0.0, 1.0] for all targets.

        Args:
            drug_concentration (float): Concentration of the drug

        Returns:
            Union[float, dict]: Inhibition percentage(s) [0.0, 1.0]
        """
        results = {name: self.calculate_occupancy(drug_concentration, name) for name in self.targets}
        if len(results) == 1 and "default" in results:
            return results["default"]
        return results

