"""
CellTypist — Fast logistic regression-based cell type annotation

Reference: https://www.celltypist.org/
Dominguez Conde, C. et al. Science (2022)
"""

from typing import Optional

import anndata
import numpy as np
import pandas as pd
import scanpy as sc


class CelltypistAnnotator:
    """CellTypist-based fast cell type prediction.

    Uses pre-trained logistic regression models for rapid annotation.
    Suitable for large-scale initial screening.

    Parameters
    ----------
    model_name : str
        CellTypist model to use, e.g., 'Immune_All_High.pkl'.
    majority_voting : bool
        Whether to use majority voting.
    """

    def __init__(
        self,
        model_name: str = "Allen_WMB.pkl",
        majority_voting: bool = True,
    ):
        self.model_name = model_name
        self.majority_voting = majority_voting

    def annotate(
        self,
        adata: anndata.AnnData,
        level: str = "subclass",
    ) -> pd.DataFrame:
        """Annotate cells using CellTypist.

        Parameters
        ----------
        adata : AnnData
            Log-normalized expression data.
        level : str
            Taxonomy level for annotation output.

        Returns
        -------
        pd.DataFrame
            Cells with CellTypist predictions and confidence scores.
        """
        n_cells = adata.n_obs
        print("logistic regression prediction...")

        try:
            import celltypist

            # Normalize to 10,000 counts per cell
            adata_norm = adata.copy()
            sc.pp.normalize_total(adata_norm, target_sum=1e4)
            sc.pp.log1p(adata_norm)

            # CellTypist needs gene symbols and count-normalized data
            # For Phase 2: use Allen WMB model
            # predictions = celltypist.annotate(
            #     adata_norm,
            #     model=self.model_name,
            #     majority_voting=self.majority_voting,
            # )

            # Phase 1: check model availability
            models_available = celltypist.models.models_description()
            print(f"    CellTypist models available: {len(models_available)}")

            result = pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "celltypist_label": ["pending_model"] * n_cells,
                "celltypist_score": [0.0] * n_cells,
            })

            print(f"    P1 framework ready; Allen model needed for Phase 2")
            return result

        except ImportError:
            print("    celltypist not installed; skipping")
            return pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "celltypist_label": ["unavailable"] * n_cells,
                "celltypist_score": [0.0] * n_cells,
            })
        except Exception as e:
            print(f"    CellTypist failed: {e}")
            return pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "celltypist_label": ["error"] * n_cells,
                "celltypist_score": [0.0] * n_cells,
            })
