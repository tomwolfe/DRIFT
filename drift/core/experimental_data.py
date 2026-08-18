"""
Experimental Dataset Registry for DRIFT Validation.

Provides structured access to published experimental datasets for
validating multi-scale drug response simulations.
"""

import os
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Union
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Default assets directory
_DEFAULT_ASSETS_DIR = Path(__file__).parent.parent.parent / "assets" / "validation"


@dataclass
class DatasetMetadata:
    """Metadata for an experimental dataset."""
    name: str
    description: str
    source: str
    doi: Optional[str]
    cell_line: str
    organism: str
    drug: str
    drug_target: str
    assay_type: str
    time_unit: str
    concentration_unit: Optional[str]
    readouts: List[str]
    file_path: Optional[str] = None
    notes: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "source": self.source,
            "doi": self.doi,
            "cell_line": self.cell_line,
            "organism": self.organism,
            "drug": self.drug,
            "drug_target": self.drug_target,
            "assay_type": self.assay_type,
            "time_unit": self.time_unit,
            "concentration_unit": self.concentration_unit,
            "readouts": self.readouts,
            "notes": self.notes,
        }


@dataclass
class ValidationDataset:
    """A loaded experimental dataset ready for comparison."""
    metadata: DatasetMetadata
    data: pd.DataFrame
    time_col: str = "time"
    readout_columns: Dict[str, str] = field(default_factory=dict)

    @property
    def time_values(self) -> np.ndarray:
        return self.data[self.time_col].values

    def get_readout(self, name: str) -> np.ndarray:
        col = self.readout_columns.get(name, name)
        return self.data[col].values

    def normalize_to_baseline(self, columns: Optional[List[str]] = None) -> "ValidationDataset":
        """Normalize readouts to their t=0 values."""
        cols = columns or [c for c in self.data.columns if c != self.time_col]
        normalized = self.data.copy()
        for col in cols:
            baseline = normalized[col].iloc[0]
            if baseline != 0:
                normalized[col] = normalized[col] / baseline
        return ValidationDataset(
            metadata=self.metadata,
            data=normalized,
            time_col=self.time_col,
            readout_columns=self.readout_columns.copy(),
        )


class DatasetRegistry:
    """
    Registry of published experimental datasets for DRIFT validation.

    Provides loaders, metadata, and unit conversion utilities for
    comparing simulation results against real experimental data.
    """

    def __init__(self, assets_dir: Optional[Union[str, Path]] = None):
        self.assets_dir = Path(assets_dir) if assets_dir else _DEFAULT_ASSETS_DIR
        self._datasets: Dict[str, DatasetMetadata] = {}
        self._register_built_in_datasets()

    def _register_built_in_datasets(self) -> None:
        """Register the built-in validation datasets."""
        self._datasets["alpelisib_mcf7"] = DatasetMetadata(
            name="alpelisib_mcf7",
            description=(
                "Alpelisib (PI3Kα inhibitor) dose-response and time-course "
                "in MCF7 breast cancer cells. Reports pAKT inhibition and "
                "cell viability over 72h."
            ),
            source="Based on Fritsch et al. 2014, Molecular Cancer Therapeutics",
            doi="10.1158/1535-7163.MCT-13-0753",
            cell_line="MCF7",
            organism="Homo sapiens",
            drug="Alpelisib (BYL719)",
            drug_target="PI3Kα",
            assay_type="Western blot (pAKT), CellTiter-Glo (viability)",
            time_unit="hours",
            concentration_unit="nM",
            readouts=["pAKT", "viability"],
            file_path=str(self.assets_dir / "alpelisib_mcf7_dose_response.csv"),
            notes="IC50 ~5 nM for pAKT inhibition, ~50 nM for viability at 72h",
        )

        self._datasets["rapamycin_time_course"] = DatasetMetadata(
            name="rapamycin_time_course",
            description=(
                "Rapamycin (mTORC1 inhibitor) time-course phosphoproteomics. "
                "Reports pS6K1, p4EBP1, and pAKT dynamics over 2 hours in "
                "HEK293 cells."
            ),
            source="Based on mTORC1 substrate kinetics (Sarbassov et al. 2004)",
            doi="10.1126/science.1097877",
            cell_line="HEK293",
            organism="Homo sapiens",
            drug="Rapamycin",
            drug_target="mTORC1",
            assay_type="Phosphoproteomics (Western blot densitometry)",
            time_unit="minutes",
            concentration_unit="nM",
            readouts=["pS6K1", "p4EBP1", "pAKT"],
            file_path=str(self.assets_dir / "rapamycin_time_course.csv"),
            notes="pS6K1 inhibited within 5 min; p4EBP1 partially resistant",
        )

        self._datasets["erlotinib_metabolic"] = DatasetMetadata(
            name="erlotinib_metabolic",
            description=(
                "Erlotinib (EGFR inhibitor) metabolic reprogramming in "
                "HCC827 NSCLC cells. Reports glucose uptake, lactate "
                "production, and oxygen consumption over 24h."
            ),
            source="Based on EGFR-TKI metabolic reprogramming data (Goggin et al. 2009)",
            doi="10.1158/1535-7163.MCT-09-0535",
            cell_line="HCC827",
            organism="Homo sapiens",
            drug="Erlotinib",
            drug_target="EGFR",
            assay_type="Seahorse XF (ECAR, OCR), enzymatic assays",
            time_unit="hours",
            concentration_unit="uM",
            readouts=["glucose_uptake", "lactate_production", "oxygen_consumption"],
            file_path=str(self.assets_dir / "erlotinib_metabolic_shift.csv"),
            notes="EGFR-TKI shifts cells from glycolysis to oxidative phosphorylation",
        )

    def list_datasets(self) -> List[str]:
        """Return names of all registered datasets."""
        return list(self._datasets.keys())

    def get_metadata(self, name: str) -> DatasetMetadata:
        """Get metadata for a registered dataset."""
        if name not in self._datasets:
            available = ", ".join(self._datasets.keys())
            raise KeyError(f"Dataset '{name}' not found. Available: {available}")
        return self._datasets[name]

    def load(self, name: str, normalize: bool = False) -> ValidationDataset:
        """
        Load a registered dataset.

        Args:
            name: Dataset name
            normalize: If True, normalize readouts to t=0 baseline

        Returns:
            ValidationDataset with loaded data
        """
        meta = self.get_metadata(name)
        if meta.file_path is None:
            raise FileNotFoundError(f"No file path registered for dataset '{name}'")

        file_path = Path(meta.file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"Data file not found: {file_path}")

        if file_path.suffix == ".csv":
            df = pd.read_csv(file_path)
        elif file_path.suffix in (".xlsx", ".xls"):
            df = pd.read_excel(file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_path.suffix}")

        dataset = ValidationDataset(
            metadata=meta,
            data=df,
            time_col="time",
            readout_columns={r: r for r in meta.readouts if r in df.columns},
        )

        if normalize:
            dataset = dataset.normalize_to_baseline()

        logger.info(f"Loaded dataset '{name}': {len(df)} rows, readouts: {list(dataset.readout_columns.keys())}")
        return dataset

    def load_from_file(
        self,
        file_path: Union[str, Path],
        metadata: Optional[DatasetMetadata] = None,
        time_col: str = "time",
        normalize: bool = False,
    ) -> ValidationDataset:
        """
        Load a dataset from an arbitrary file with optional metadata.

        Args:
            file_path: Path to CSV/Excel file
            metadata: Optional metadata (auto-detected from filename if None)
            time_col: Name of the time column
            normalize: If True, normalize readouts to t=0 baseline

        Returns:
            ValidationDataset
        """
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        if file_path.suffix == ".csv":
            df = pd.read_csv(file_path)
        elif file_path.suffix in (".xlsx", ".xls"):
            df = pd.read_excel(file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_path.suffix}")

        if time_col not in df.columns:
            raise ValueError(f"Time column '{time_col}' not found. Available: {list(df.columns)}")

        readout_columns = {c: c for c in df.columns if c != time_col}

        if metadata is None:
            metadata = DatasetMetadata(
                name=file_path.stem,
                description=f"User-provided dataset from {file_path.name}",
                source="User-provided",
                doi=None,
                cell_line="unknown",
                organism="unknown",
                drug="unknown",
                drug_target="unknown",
                assay_type="unknown",
                time_unit="unknown",
                concentration_unit=None,
                readouts=list(readout_columns.keys()),
                file_path=str(file_path),
            )

        dataset = ValidationDataset(
            metadata=metadata,
            data=df,
            time_col=time_col,
            readout_columns=readout_columns,
        )

        if normalize:
            dataset = dataset.normalize_to_baseline()

        return dataset

    def register_dataset(self, metadata: DatasetMetadata) -> None:
        """Register a new dataset metadata entry."""
        self._datasets[metadata.name] = metadata
        logger.info(f"Registered dataset: {metadata.name}")

    @staticmethod
    def convert_time(values: np.ndarray, from_unit: str, to_unit: str) -> np.ndarray:
        """
        Convert time values between units.

        Supported units: seconds, minutes, hours, days
        """
        to_seconds = {
            "seconds": 1.0,
            "sec": 1.0,
            "s": 1.0,
            "minutes": 60.0,
            "min": 60.0,
            "m": 60.0,
            "hours": 3600.0,
            "hr": 3600.0,
            "h": 3600.0,
            "days": 86400.0,
            "d": 86400.0,
        }

        from_key = from_unit.lower()
        to_key = to_unit.lower()

        if from_key not in to_seconds:
            raise ValueError(f"Unknown time unit: {from_unit}")
        if to_key not in to_seconds:
            raise ValueError(f"Unknown time unit: {to_unit}")

        return values * (to_seconds[from_key] / to_seconds[to_key])

    @staticmethod
    def convert_concentration(values: np.ndarray, from_unit: str, to_unit: str) -> np.ndarray:
        """
        Convert concentration values between units.

        Supported units: M, mM, uM, nM, pM
        """
        to_molar = {
            "M": 1.0,
            "mM": 1e-3,
            "um": 1e-6,
            "uM": 1e-6,
            "nm": 1e-9,
            "nM": 1e-9,
            "pm": 1e-12,
            "pM": 1e-12,
        }

        from_key = from_unit.lower()
        to_key = to_unit.lower()

        if from_key not in to_molar:
            raise ValueError(f"Unknown concentration unit: {from_unit}")
        if to_key not in to_molar:
            raise ValueError(f"Unknown concentration unit: {to_unit}")

        return values * (to_molar[from_key] / to_molar[to_key])
