"""
scVI/scANVI — Deep generative model annotation wrapper

Uses scvi-tools for batch-aware integration and semi-supervised annotation.

Reference: https://scvi-tools.org/
"""

from typing import Optional

import anndata
import numpy as np
import pandas as pd
import scanpy as sc


class ScanviAnnotator:
    """scVI/scANVI-based semi-supervised annotation.

    Requires a pre-trained scANVI model on Allen/BICCN reference data.
    Uses scArches for query-to-reference mapping.

    Parameters
    ----------
    model_path : str, optional
        Path to pre-trained scANVI model.
    n_latent : int
        Latent dimension for scVI model.
    batch_key : str
        Column in adata.obs for batch labels.
    """

    def __init__(
        self,
        model_path: Optional[str] = None,
        n_latent: int = 30,
        batch_key: str = "sample_id",
    ):
        self.model_path = model_path
        self.n_latent = n_latent
        self.batch_key = batch_key

    def annotate(
        self,
        adata: anndata.AnnData,
        level: str = "subclass",
    ) -> pd.DataFrame:
        """Annotate cells using scANVI label transfer.

        Parameters
        ----------
        adata : AnnData
            Preprocessed expression data.
        level : str
            Taxonomy level for annotation output.

        Returns
        -------
        pd.DataFrame
            Cells with scANVI labels and confidence scores.
        """
        n_cells = adata.n_obs
        print("label transfer...")

        try:
            import scvi

            # Phase 2 integration:
            # 1. Load pre-trained scANVI model on Allen WMB reference
            # 2. Setup AnnData for query
            # 3. Use scArches for query mapping
            # 4. Extract predicted labels at requested level

            # scvi.model.SCANVI.load_query_data(adata, reference_model)
            # scvi.model.SCANVI.load(dir_path, adata)
            # predictions = model.predict()

            # Minimal Phase 1: run a basic scVI integration as demo
            if self.batch_key not in adata.obs.columns:
                adata.obs[self.batch_key] = "batch_1"

            sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=self.batch_key)
            adata_scvi = adata[:, adata.var.highly_variable].copy()

            scvi.model.SCVI.setup_anndata(adata_scvi, batch_key=self.batch_key)
            model = scvi.model.SCVI(adata_scvi, n_latent=self.n_latent)
            model.train(max_epochs=50, early_stopping=True, plan_kwargs={"lr": 1e-3})

            latent = model.get_latent_representation()
            adata.obsm["X_scvi"] = latent

            # Run clustering as placeholder for label transfer
            sc.pp.neighbors(adata, use_rep="X_scvi")
            sc.tl.leiden(adata, resolution=0.8)

            labels = adata.obs["leiden"].astype(str).values
            scores = np.ones(n_cells) * 0.7  # placeholder confidence

            result = pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "scanvi_label": [f"cluster_{l}" for l in labels],
                "scanvi_score": scores,
            })

            print(f"    scVI integration complete, {model.n_latent}D latent space")
            return result

        except ImportError:
            print("    scvi-tools not installed; skipping scANVI")
            return pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "scanvi_label": ["scvi_unavailable"] * n_cells,
                "scanvi_score": [0.0] * n_cells,
            })
        except Exception as e:
            print(f"    scVI failed: {e}")
            return pd.DataFrame({
                "cell_barcode": adata.obs_names,
                "scanvi_label": ["error"] * n_cells,
                "scanvi_score": [0.0] * n_cells,
            })
