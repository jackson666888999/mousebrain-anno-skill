"""
SingleR — Reference-based correlation annotation

Reference: https://www.bioconductor.org/packages/SingleR/
Aran, D. et al. Nature Immunology (2019)
"""

from typing import Optional

import anndata
import numpy as np
import pandas as pd


class SinglerAnnotator:
    """SingleR-based correlation annotation.

    Uses reference expression profiles for cell type identification
    via Spearman correlation.

    Parameters
    ----------
    ref_data : str, optional
        Path to reference dataset (h5ad).
    ref_labels : str, optional
        Path to reference labels (CSV).
    """

    def __init__(
        self,
        ref_data: Optional[str] = None,
        ref_labels: Optional[str] = None,
    ):
        self.ref_data = ref_data
        self.ref_labels = ref_labels

    def annotate(
        self,
        adata: anndata.AnnData,
        level: str = "subclass",
    ) -> pd.DataFrame:
        """Annotate cells using SingleR correlation.

        Parameters
        ----------
        adata : AnnData
            Log-normalized expression data.
        level : str
            Taxonomy level for annotation output.

        Returns
        -------
        pd.DataFrame
            Cells with SingleR predictions and scores.
        """
        n_cells = adata.n_obs
        print("correlation-based annotation...")

        # SingleR originally runs in R (Bioconductor).
        # We provide:
        # 1. R bridge via rpy2 (if available)
        # 2. Python-native correlation-based classifier (Phase 2)

        try:
            # Attempt rpy2-based R interop
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri

            pandas2ri.activate()

            print(f"    R interop via rpy2 available")
            # Full R pipeline:
            # ro.r('library(SingleR)')
            # ro.r('library(celldex)')
            # ref = ro.r('MouseRNAseqData()')
            # pred = ro.r('SingleR(test=query, ref=ref, labels=ref$label.main)')

            result = pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "singler_label": ["r_interop_pending"] * n_cells,
                "singler_score": [0.0] * n_cells,
            })

            return result

        except ImportError:
            # Python-native fallback: simple correlation
            print(f"    rpy2 unavailable; Python-native mode not implemented in P1")
            return pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "singler_label": ["unavailable"] * n_cells,
                "singler_score": [0.0] * n_cells,
            })
        except Exception as e:
            print(f"    SingleR failed: {e}")
            return pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "singler_label": ["error"] * n_cells,
                "singler_score": [0.0] * n_cells,
            })
