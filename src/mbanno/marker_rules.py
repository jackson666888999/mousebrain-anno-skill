"""
Marker Rules — Rule-based cell type annotation using known mouse brain markers

Sources: Allen WMB Atlas (Nature 2023), published literature.
Markers are derived from public data, not copied from paper tables.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata
import numpy as np
import pandas as pd


class MarkerRuleEngine:
    """Rule-based annotation using marker gene combinations.

    Provides interpretable, non-black-box cell type calls based on
    well-established mouse brain marker genes.

    Parameters
    ----------
    marker_dir : str or Path
        Directory containing marker gene definition files.
    """

    def __init__(self, marker_dir: str = "data/references/markers"):
        self.marker_dir = Path(marker_dir)
        self.rules = self._load_rules()

    def _load_rules(self) -> Dict:
        """Load marker gene rules from disk or built-in defaults."""
        marker_file = self.marker_dir / "marker_genes.json"
        if marker_file.exists():
            with open(marker_file) as f:
                return json.load(f)
        return self._get_default_rules()

    def _get_default_rules(self) -> Dict:
        """Built-in default mouse brain marker rules."""
        return {
            "neuron": {
                "glutamatergic": ["Slc17a7", "Slc17a6", "Camk2a", "Neurod6"],
                "GABAergic": ["Gad1", "Gad2", "Slc32a1"],
                "cortical_IT": ["Satb2", "Cux1", "Cux2"],
                "cortical_ET": ["Fezf2", "Bcl11b", "Foxp2"],
                "cortical_CT": ["Tle4", "Foxp2"],
                "MSN_D1": ["Drd1", "Tac1", "Ppp1r1b"],
                "MSN_D2": ["Drd2", "Adora2a", "Penk"],
                "cholinergic": ["Chat", "Slc5a7"],
                "dopaminergic": ["Th", "Slc6a3", "Ddc"],
                "serotonergic": ["Tph2", "Slc6a4"],
            },
            "glia": {
                "astrocyte": ["Aqp4", "Aldh1l1", "Slc1a3", "Gfap"],
                "OPC": ["Pdgfra", "Cspg4", "Sox10"],
                "oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mag"],
                "microglia": ["P2ry12", "Tmem119", "Cx3cr1", "Csf1r"],
            },
            "vascular": {
                "endothelial": ["Cldn5", "Pecam1", "Flt1", "Kdr"],
                "pericyte": ["Pdgfrb", "Rgs5"],
                "smooth_muscle": ["Acta2", "Tagln"],
            },
            "other": {
                "ependymal": ["Foxj1"],
                "choroid_plexus": ["Ttr", "Krt8", "Krt18"],
            },
            "data_source": "Allen WMB Atlas (Nature 2023) + literature",
            "license": "CC-BY-4.0 (Allen data); Public domain (literature)",
        }

    def annotate(
        self,
        adata: anndata.AnnData,
    ) -> pd.DataFrame:
        """Score and classify cells using marker gene rules.

        Parameters
        ----------
        adata : AnnData
            Log-normalized expression data.

        Returns
        -------
        pd.DataFrame
            Per-cell marker scores and predicted labels.
        """
        results = []

        for category, cell_types in self.rules.items():
            if category in ("data_source", "license"):
                continue

            for cell_type, markers in cell_types.items():
                valid_markers = [g for g in markers if g in adata.var_names]
                if not valid_markers:
                    continue

                # Mean expression of valid markers per cell
                marker_idx = adata.var_names.isin(valid_markers)
                expr = adata.X[:, marker_idx]
                if hasattr(expr, "toarray"):
                    expr = expr.toarray()

                mean_expr = np.mean(expr, axis=1)

                results.append({
                    "category": category,
                    "cell_type": cell_type,
                    "markers_used": len(valid_markers),
                    "markers_total": len(markers),
                    "mean_expression": float(np.mean(mean_expr)),
                    "expression_per_cell": mean_expr,
                })

        df = pd.DataFrame(results)
        return df

    def classify_cells(
        self,
        adata: anndata.AnnData,
        min_percentile: float = 0.5,
    ) -> pd.DataFrame:
        """Assign marker-based labels to individual cells.

        Parameters
        ----------
        adata : AnnData
            Log-normalized expression data.
        min_percentile : float
            Minimum expression percentile to call a type.

        Returns
        -------
        pd.DataFrame
            Cell barcodes, best marker label, and scores.
        """
        all_scores = {}

        for category, cell_types in self.rules.items():
            if category in ("data_source", "license"):
                continue

            for cell_type, markers in cell_types.items():
                valid_markers = [g for g in markers if g in adata.var_names]
                if not valid_markers:
                    continue

                marker_idx = adata.var_names.isin(valid_markers)
                expr = adata.X[:, marker_idx]
                if hasattr(expr, "toarray"):
                    expr = expr.toarray()

                scores = np.mean(expr, axis=1)
                all_scores[cell_type] = scores

        if not all_scores:
            return pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "marker_label": ["unassigned"] * adata.n_obs,
                "marker_score": [0.0] * adata.n_obs,
            })

        # Stack all scores into a matrix [n_cells x n_types]
        types = list(all_scores.keys())
        score_matrix = np.column_stack([all_scores[t] for t in types])

        # Best type per cell
        best_idx = np.argmax(score_matrix, axis=1)
        best_scores = score_matrix[np.arange(adata.n_obs), best_idx]

        # Compute margin (best - second best)
        sorted_indices = np.argsort(-score_matrix, axis=1)
        if score_matrix.shape[1] >= 2:
            second_best = score_matrix[np.arange(adata.n_obs), sorted_indices[:, 1]]
            margins = best_scores - second_best
        else:
            margins = best_scores

        result = pd.DataFrame({
            "cell_barcode": adata.obs_names,
            "marker_label": [types[i] for i in best_idx],
            "marker_score": best_scores,
            "marker_margin": margins,
        })

        return result
