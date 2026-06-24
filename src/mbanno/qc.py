"""
QC — Quality control and preprocessing for mouse brain single-cell data
"""

from typing import Optional, Dict

import anndata
import numpy as np
import scanpy as sc


class QCProcessor:
    """Standardized QC for mouse brain single-cell/nucleus RNA-seq data.

    Parameters
    ----------
    min_genes : int
        Minimum number of genes detected per cell.
    min_counts : int
        Minimum total UMI counts per cell.
    max_pct_mt : float
        Maximum mitochondrial gene percentage.
    doublet_threshold : float
        Scrublet doublet score threshold.
    """

    def __init__(
        self,
        min_genes: int = 200,
        min_counts: int = 500,
        max_pct_mt: float = 20.0,
        doublet_threshold: float = 0.25,
    ):
        self.min_genes = min_genes
        self.min_counts = min_counts
        self.max_pct_mt = max_pct_mt
        self.doublet_threshold = doublet_threshold

    def run(self, adata: anndata.AnnData, detect_doublets: bool = True) -> anndata.AnnData:
        """Run full QC pipeline.

        Parameters
        ----------
        adata : AnnData
            Raw count matrix.
        detect_doublets : bool
            Whether to run Scrublet doublet detection.

        Returns
        -------
        AnnData
            QC-filtered data.
        """
        print(f"  [QC] Input cells: {adata.n_obs:,}")

        # Mitochondrial genes (mouse naming: mt-)
        adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")

        # Ribosomal genes
        adata.var["ribo"] = adata.var_names.str.lower().str.startswith(("rps", "rpl"))

        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
        )

        # Doublet detection via Scrublet
        if detect_doublets:
            try:
                import scrublet as scr
                scrub = scr.Scrublet(adata.X)
                doublet_scores, predicted_doublets = scrub.scrub_doublets()
                adata.obs["doublet_score"] = doublet_scores
                adata.obs["predicted_doublet"] = predicted_doublets
                print(f"  [QC] Doublet rate (est.): {predicted_doublets.mean():.2%}")
            except ImportError:
                print("  [QC] scrublet not installed; skipping doublet detection.")
                adata.obs["doublet_score"] = 0.0
                adata.obs["predicted_doublet"] = False

        # Apply filters
        n_before = adata.n_obs

        keep = (
            (adata.obs["n_genes_by_counts"] >= self.min_genes)
            & (adata.obs["total_counts"] >= self.min_counts)
            & (adata.obs["pct_counts_mt"] < self.max_pct_mt)
        )

        if "doublet_score" in adata.obs:
            keep = keep & (adata.obs["doublet_score"] < self.doublet_threshold)

        adata = adata[keep, :].copy()

        n_after = adata.n_obs
        print(f"  [QC] After filtering: {n_after:,} cells "
              f"(removed {n_before - n_after:,})")

        # Normalization
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Highly variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
        n_hvg = adata.var["highly_variable"].sum()
        print(f"  [QC] HVGs: {n_hvg}")

        return adata

    def get_summary(self, adata: anndata.AnnData) -> Dict:
        """Generate QC summary statistics."""
        return {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "median_genes_per_cell": float(np.median(adata.obs["n_genes_by_counts"])),
            "median_counts_per_cell": float(np.median(adata.obs["total_counts"])),
            "median_pct_mt": float(np.median(adata.obs["pct_counts_mt"])),
            "median_pct_ribo": float(np.median(adata.obs.get("pct_counts_ribo", [0]))),
            "doublet_rate": float(adata.obs.get("predicted_doublet", [False]).mean()),
        }
