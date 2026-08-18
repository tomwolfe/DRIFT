"""
Validation module for DRIFT.

Provides tools for comparing simulation results against experimental data,
including multiple statistical metrics and structured reporting.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Union, Tuple
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


@dataclass
class MetricResult:
    """Result for a single comparison metric."""
    name: str
    value: float
    passed: bool
    threshold: Optional[float] = None
    details: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "value": self.value,
            "passed": self.passed,
            "threshold": self.threshold,
            "details": self.details,
        }


@dataclass
class SpeciesComparison:
    """Comparison results for a single species/readout."""
    species_name: str
    experiment_values: np.ndarray
    simulation_values: np.ndarray
    time_values: np.ndarray
    metrics: Dict[str, MetricResult] = field(default_factory=dict)

    @property
    def overall_passed(self) -> bool:
        return all(m.passed for m in self.metrics.values())


@dataclass
class ValidationReport:
    """Structured report from a validation run."""
    dataset_name: str
    comparisons: Dict[str, SpeciesComparison] = field(default_factory=dict)
    overall_passed: bool = False
    summary: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "dataset_name": self.dataset_name,
            "overall_passed": self.overall_passed,
            "summary": self.summary,
            "metadata": self.metadata,
            "comparisons": {
                name: {
                    "species": comp.species_name,
                    "passed": comp.overall_passed,
                    "metrics": {m_name: m.to_dict() for m_name, m in comp.metrics.items()},
                }
                for name, comp in self.comparisons.items()
            },
        }


def compute_rmse(sim: np.ndarray, exp: np.ndarray) -> float:
    """Root Mean Square Error."""
    return float(np.sqrt(np.mean((sim - exp) ** 2)))


def compute_r_squared(sim: np.ndarray, exp: np.ndarray) -> float:
    """Coefficient of determination (R²)."""
    ss_res = np.sum((exp - sim) ** 2)
    ss_tot = np.sum((exp - np.mean(exp)) ** 2)
    if ss_tot == 0:
        return 0.0
    return float(1.0 - (ss_res / ss_tot))


def compute_spearman(sim: np.ndarray, exp: np.ndarray) -> float:
    """Spearman rank correlation coefficient."""
    from scipy.stats import spearmanr
    corr, _ = spearmanr(sim, exp)
    return float(corr) if np.isfinite(corr) else 0.0


def compute_bland_altman(sim: np.ndarray, exp: np.ndarray) -> Dict[str, float]:
    """
    Bland-Altman limits of agreement.

    Returns dict with 'mean_diff', 'std_diff', 'lower_limit', 'upper_limit'.
    """
    diffs = sim - exp
    mean_diff = float(np.mean(diffs))
    std_diff = float(np.std(diffs, ddof=1))
    return {
        "mean_diff": mean_diff,
        "std_diff": std_diff,
        "lower_limit": mean_diff - 1.96 * std_diff,
        "upper_limit": mean_diff + 1.96 * std_diff,
    }


def compute_fold_change_agreement(sim: np.ndarray, exp: np.ndarray) -> float:
    """
    Fraction of time points where simulation and experiment agree on
    direction of change (fold-change sign agreement).
    """
    if len(sim) < 2:
        return 1.0
    sim_direction = np.sign(np.diff(sim))
    exp_direction = np.sign(np.diff(exp))
    agreement = np.mean(sim_direction == exp_direction)
    return float(agreement)


class ExternalValidator:
    """
    Framework for importing and benchmarking simulation results against
    external experimental datasets.

    Supports both file-based and DatasetRegistry-based validation with
    comprehensive statistical metrics.
    """

    def __init__(self, experimental_data_path: Optional[str] = None):
        """
        Initialize the validator.

        Args:
            experimental_data_path: Path to a CSV/Excel file containing experimental data.
        """
        self.exp_data = None
        if experimental_data_path:
            self.load_experimental_data(experimental_data_path)

    def load_experimental_data(self, path: str) -> None:
        """Loads experimental data from a file."""
        if path.endswith(".csv"):
            self.exp_data = pd.read_csv(path)
        elif path.endswith(".xlsx") or path.endswith(".xls"):
            self.exp_data = pd.read_excel(path)
        else:
            raise ValueError("Unsupported file format. Use CSV or Excel.")
        logger.info(f"Loaded experimental data from {path}")

    def load_dataset(self, dataset: Any) -> None:
        """
        Load data from a ValidationDataset object.

        Args:
            dataset: A ValidationDataset from DatasetRegistry
        """
        self.exp_data = dataset.data.copy()
        self._time_col = dataset.time_col
        logger.info(f"Loaded dataset '{dataset.metadata.name}': {len(self.exp_data)} rows")

    def benchmark(
        self,
        simulation_history: Dict[str, Any],
        mapping: Dict[str, str],
        time_col: str = "time",
    ) -> Dict[str, Any]:
        """
        Compare simulation history with experimental data.

        Args:
            simulation_history: The history dict from a simulation run.
            mapping: Mapping from experimental column names to simulation species/keys.
                     e.g., {"pAKT_level": "AKT", "OCR": "growth"}
            time_col: The column name for time in experimental data.

        Returns:
            dict: Benchmark metrics (MSE, correlation) for each mapped entity.
        """
        if self.exp_data is None:
            raise ValueError("No experimental data loaded.")

        results: Dict[str, Any] = {}
        sim_time = simulation_history["time"]

        from scipy.interpolate import CubicSpline

        for exp_col, sim_key in mapping.items():
            if exp_col not in self.exp_data.columns:
                logger.warning(f"Column {exp_col} not found in experimental data.")
                continue

            try:
                sorted_exp = self.exp_data.sort_values(time_col)
                cs = CubicSpline(sorted_exp[time_col], sorted_exp[exp_col], extrapolate=True)
                exp_val = cs(sim_time)
                if (sorted_exp[exp_col] >= 0).all():
                    exp_val = np.maximum(0, exp_val)
            except Exception as e:
                logger.warning(f"CubicSpline interpolation failed for {exp_col}: {e}. Falling back to linear.")
                exp_val = np.interp(
                    sim_time,
                    self.exp_data[time_col],
                    self.exp_data[exp_col],
                )

            if sim_key in simulation_history:
                sim_val = simulation_history[sim_key]
            elif "signaling" in simulation_history:
                species = simulation_history.get("species_names", ["PI3K", "AKT", "mTOR"])
                if sim_key in species:
                    idx = species.index(sim_key)
                    sim_val = simulation_history["signaling"][:, idx]
                else:
                    logger.warning(f"Species {sim_key} not found in signaling history.")
                    continue
            else:
                logger.warning(f"Key {sim_key} not found in simulation history.")
                continue

            mse = float(np.mean((sim_val - exp_val) ** 2))
            corr = float(np.corrcoef(sim_val, exp_val)[0, 1])

            results[sim_key] = {
                "mse": mse,
                "correlation": corr,
            }

        return results

    def full_benchmark(
        self,
        simulation_history: Dict[str, Any],
        mapping: Dict[str, str],
        time_col: str = "time",
        species_names: Optional[List[str]] = None,
        pass_thresholds: Optional[Dict[str, float]] = None,
    ) -> ValidationReport:
        """
        Comprehensive comparison with structured report.

        Args:
            simulation_history: The history dict from a simulation run.
            mapping: Mapping from experimental column names to simulation species/keys.
            time_col: Column name for time in experimental data.
            species_names: List of species names for signaling extraction.
            pass_thresholds: Dict of metric_name -> threshold for pass/fail.

        Returns:
            ValidationReport with all metrics
        """
        if self.exp_data is None:
            raise ValueError("No experimental data loaded.")

        if pass_thresholds is None:
            pass_thresholds = {
                "rmse": 0.15,
                "r_squared": 0.7,
                "spearman": 0.6,
                "direction_agreement": 0.6,
            }

        sim_time = simulation_history["time"]
        if species_names is None:
            species_names = simulation_history.get("species_names", ["PI3K", "AKT", "mTOR"])

        report = ValidationReport(
            dataset_name="user_dataset",
            metadata={
                "n_timepoints": len(sim_time),
                "time_range": (float(sim_time[0]), float(sim_time[-1])),
                "mapping": mapping,
            },
        )

        from scipy.interpolate import CubicSpline

        for exp_col, sim_key in mapping.items():
            if exp_col not in self.exp_data.columns:
                logger.warning(f"Column {exp_col} not found in experimental data.")
                continue

            try:
                sorted_exp = self.exp_data.sort_values(time_col)
                cs = CubicSpline(sorted_exp[time_col], sorted_exp[exp_col], extrapolate=True)
                exp_val = cs(sim_time)
                if (sorted_exp[exp_col] >= 0).all():
                    exp_val = np.maximum(0, exp_val)
            except Exception as e:
                logger.warning(f"CubicSpline failed for {exp_col}: {e}. Using linear fallback.")
                exp_val = np.interp(sim_time, self.exp_data[time_col], self.exp_data[exp_col])

            if sim_key in simulation_history:
                sim_val = simulation_history[sim_key]
            elif "signaling" in simulation_history and sim_key in species_names:
                idx = species_names.index(sim_key)
                sim_val = simulation_history["signaling"][:, idx]
            else:
                logger.warning(f"Key '{sim_key}' not found in simulation history.")
                continue

            # Compute all metrics
            rmse = compute_rmse(sim_val, exp_val)
            r2 = compute_r_squared(sim_val, exp_val)
            spearman = compute_spearman(sim_val, exp_val)
            bland_altman = compute_bland_altman(sim_val, exp_val)
            direction = compute_fold_change_agreement(sim_val, exp_val)

            metrics = {
                "rmse": MetricResult(
                    name="RMSE",
                    value=rmse,
                    passed=rmse <= pass_thresholds.get("rmse", 0.15),
                    threshold=pass_thresholds.get("rmse", 0.15),
                    details=f"Root mean square error: {rmse:.4f}",
                ),
                "r_squared": MetricResult(
                    name="R²",
                    value=r2,
                    passed=r2 >= pass_thresholds.get("r_squared", 0.7),
                    threshold=pass_thresholds.get("r_squared", 0.7),
                    details=f"Coefficient of determination: {r2:.4f}",
                ),
                "spearman": MetricResult(
                    name="Spearman",
                    value=spearman,
                    passed=abs(spearman) >= pass_thresholds.get("spearman", 0.6),
                    threshold=pass_thresholds.get("spearman", 0.6),
                    details=f"Spearman rank correlation: {spearman:.4f}",
                ),
                "direction_agreement": MetricResult(
                    name="Direction Agreement",
                    value=direction,
                    passed=direction >= pass_thresholds.get("direction_agreement", 0.6),
                    threshold=pass_thresholds.get("direction_agreement", 0.6),
                    details=f"Fraction of time points with correct direction: {direction:.2%}",
                ),
                "bland_altman_mean": MetricResult(
                    name="Bland-Altman Mean Diff",
                    value=bland_altman["mean_diff"],
                    passed=abs(bland_altman["mean_diff"]) < pass_thresholds.get("rmse", 0.15),
                    details=(
                        f"Mean difference: {bland_altman['mean_diff']:.4f}, "
                        f"LoA: [{bland_altman['lower_limit']:.4f}, {bland_altman['upper_limit']:.4f}]"
                    ),
                ),
            }

            comp = SpeciesComparison(
                species_name=sim_key,
                experiment_values=exp_val,
                simulation_values=sim_val,
                time_values=sim_time,
                metrics=metrics,
            )
            report.comparisons[exp_col] = comp

        report.overall_passed = all(c.overall_passed for c in report.comparisons.values()) if report.comparisons else False
        n_passed = sum(1 for c in report.comparisons.values() if c.overall_passed)
        report.summary = f"{n_passed}/{len(report.comparisons)} readouts passed all thresholds"

        return report

    def plot_comparison(
        self,
        simulation_history: Dict[str, Any],
        mapping: Dict[str, str],
        time_col: str = "time",
        output_path: Optional[str] = None,
        report: Optional[ValidationReport] = None,
    ) -> None:
        """Visualizes the comparison between simulation and experiment."""
        if self.exp_data is None:
            raise ValueError("No experimental data loaded.")

        n_plots = len(mapping)
        fig, axes = plt.subplots(n_plots, 1, figsize=(10, 4 * n_plots), sharex=True)
        if n_plots == 1:
            axes = [axes]

        sim_time = simulation_history["time"]
        species_names = simulation_history.get("species_names", ["PI3K", "AKT", "mTOR"])

        for i, (exp_col, sim_key) in enumerate(mapping.items()):
            ax = axes[i]

            ax.scatter(self.exp_data[time_col], self.exp_data[exp_col], color="red", label="Experimental", zorder=5)

            if sim_key == "growth":
                sim_val = simulation_history["growth"]
            elif sim_key in species_names:
                idx = species_names.index(sim_key)
                sim_val = simulation_history["signaling"][:, idx]
            else:
                continue

            ax.plot(sim_time, sim_val, label="Simulation", color="blue")

            # Add metric annotations if report is provided
            if report and exp_col in report.comparisons:
                comp = report.comparisons[exp_col]
                r2 = comp.metrics.get("r_squared")
                rmse = comp.metrics.get("rmse")
                annotation = ""
                if r2:
                    annotation += f"R²={r2.value:.3f}  "
                if rmse:
                    annotation += f"RMSE={rmse.value:.4f}"
                if annotation:
                    ax.annotate(
                        annotation,
                        xy=(0.02, 0.95),
                        xycoords="axes fraction",
                        fontsize=9,
                        verticalalignment="top",
                        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
                    )

            ax.set_ylabel(sim_key)
            ax.legend()

        axes[-1].set_xlabel("Time")
        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            logger.info(f"Comparison plot saved to {output_path}")
        else:
            plt.show()

    def plot_bland_altman(
        self,
        simulation_history: Dict[str, Any],
        mapping: Dict[str, str],
        time_col: str = "time",
        output_path: Optional[str] = None,
    ) -> None:
        """Generates Bland-Altman plots for each mapped species."""
        if self.exp_data is None:
            raise ValueError("No experimental data loaded.")

        n_plots = len(mapping)
        fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 5))
        if n_plots == 1:
            axes = [axes]

        sim_time = simulation_history["time"]
        species_names = simulation_history.get("species_names", ["PI3K", "AKT", "mTOR"])

        from scipy.interpolate import CubicSpline

        for i, (exp_col, sim_key) in enumerate(mapping.items()):
            ax = axes[i]
            if exp_col not in self.exp_data.columns:
                continue

            try:
                sorted_exp = self.exp_data.sort_values(time_col)
                cs = CubicSpline(sorted_exp[time_col], sorted_exp[exp_col], extrapolate=True)
                exp_val = cs(sim_time)
                if (sorted_exp[exp_col] >= 0).all():
                    exp_val = np.maximum(0, exp_val)
            except Exception:
                exp_val = np.interp(sim_time, self.exp_data[time_col], self.exp_data[exp_col])

            if sim_key in simulation_history:
                sim_val = simulation_history[sim_key]
            elif sim_key in species_names:
                idx = species_names.index(sim_key)
                sim_val = simulation_history["signaling"][:, idx]
            else:
                continue

            mean_vals = (sim_val + exp_val) / 2.0
            diffs = sim_val - exp_val
            mean_diff = np.mean(diffs)
            std_diff = np.std(diffs, ddof=1)

            ax.scatter(mean_vals, diffs, alpha=0.6, s=30)
            ax.axhline(y=mean_diff, color="gray", linestyle="--", label=f"Mean diff: {mean_diff:.4f}")
            ax.axhline(y=mean_diff + 1.96 * std_diff, color="red", linestyle=":", label="±1.96 SD")
            ax.axhline(y=mean_diff - 1.96 * std_diff, color="red", linestyle=":")
            ax.axhline(y=0, color="black", linewidth=0.5)
            ax.set_xlabel(f"Mean of ({sim_key}, {exp_col})")
            ax.set_ylabel(f"Difference ({sim_key} - {exp_col})")
            ax.set_title(f"Bland-Altman: {sim_key}")
            ax.legend(fontsize=8)

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            logger.info(f"Bland-Altman plot saved to {output_path}")
        else:
            plt.show()
