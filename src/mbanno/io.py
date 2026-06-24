"""
I/O — Data loading, saving, and format handling
"""

import json
import warnings
from pathlib import Path
from typing import Optional, Union, Dict, Any

import anndata
import numpy as np
import pandas as pd
import yaml


def load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
    """Load a YAML config file."""
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_data_sources(manifest_path: Union[str, Path] = "data_sources.yaml") -> Dict[str, Any]:
    """Load the data sources manifest."""
    manifest_path = Path(manifest_path)
    if not manifest_path.exists():
        # Fallback: look relative to the package
        alt_path = Path(__file__).parent.parent.parent / "data_sources.yaml"
        if alt_path.exists():
            manifest_path = alt_path
        else:
            raise FileNotFoundError(
                f"Data sources manifest not found: {manifest_path}. "
                "Run `mbanno download-reference --list` to see available sources."
            )
    return load_config(manifest_path)


def read_adata(path: Union[str, Path]) -> anndata.AnnData:
    """Read an AnnData object from h5ad or h5 file."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    if path.suffix in (".h5ad", ".h5"):
        return anndata.read_h5ad(path)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}. Expected .h5ad or .h5")


def write_adata(adata: anndata.AnnData, path: Union[str, Path]):
    """Write AnnData to h5ad."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    print(f"  Saved: {path}")


def read_markers(marker_path: Union[str, Path]) -> dict:
    """Read marker gene file (JSON or TSV)."""
    marker_path = Path(marker_path)
    if marker_path.suffix == ".json":
        with open(marker_path) as f:
            return json.load(f)
    elif marker_path.suffix in (".tsv", ".csv"):
        sep = "\t" if marker_path.suffix == ".tsv" else ","
        return pd.read_csv(marker_path, sep=sep).to_dict(orient="list")
    else:
        raise ValueError(f"Unsupported marker format: {marker_path.suffix}")


def write_annotations(
    results: pd.DataFrame,
    path: Union[str, Path],
):
    """Write annotation results to TSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(path, sep="\t", index=False)
    print(f"  Saved annotations: {path}")


def write_confidence(
    confidence_df: pd.DataFrame,
    path: Union[str, Path],
):
    """Write confidence metrics to TSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    confidence_df.to_csv(path, sep="\t", index=False)
    print(f"  Saved confidence: {path}")


def write_json(data: dict, path: Union[str, Path]):
    """Write dict to JSON file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def _convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pd.Series):
            return obj.tolist()
        if isinstance(obj, pd.DataFrame):
            return obj.to_dict(orient="records")
        return obj

    cleaned = _convert(data)
    with open(path, "w") as f:
        json.dump(cleaned, f, indent=2, default=str)
