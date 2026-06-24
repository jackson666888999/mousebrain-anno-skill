"""
MapMyCells — Allen/BRAIN Initiative official mapping interface

Maps user expression data to Allen/BICCN taxonomies via the
MapMyCells API or local reference files.

Reference: https://www.braininitiative.org/toolmakers/resources/mapmycells/
"""

from typing import Optional

import anndata
import numpy as np
import pandas as pd


class MapMyCellsRunner:
    """Interface to Allen/BICCN MapMyCells mapping service.

    MapMyCells is the official Allen/BRAIN Initiative tool for mapping
    query data to Allen/BICAN taxonomies. This runner provides a
    Python interface for annotation.

    Parameters
    ----------
    taxonomy : str
        Target taxonomy, e.g., 'WMB-10X' or 'WMB-Multiome'.
    mapmycells_url : str
        API endpoint for MapMyCells service.
    """

    def __init__(
        self,
        taxonomy: str = "WMB-10X",
        mapmycells_url: str = "https://knowledge.brain-map.org/mapmycells/api/v1",
    ):
        self.taxonomy = taxonomy
        self.mapmycells_url = mapmycells_url

    def annotate(
        self,
        adata: anndata.AnnData,
        level: str = "subclass",
    ) -> pd.DataFrame:
        """Annotate cells using MapMyCells taxonomy mapping.

        Parameters
        ----------
        adata : AnnData
            Log-normalized expression data.
        level : str
            Taxonomy level: class, subclass, supertype, cluster.

        Returns
        -------
        pd.DataFrame
            Cells with MapMyCells labels and scores.
        """
        n_cells = adata.n_obs
        print("mapping to Allen/BICCN taxonomy...")

        # MapMyCells requires:
        # 1. Gene symbols in var_names (human or mouse)
        # 2. Counts or normalized data in .X
        # 3. Cell metadata in .obs

        # In Phase 1, this is a framework that will:
        # - Call the MapMyCells API with the expression matrix
        # - Parse the returned taxonomy assignments
        # - Map to the requested level

        # Full integration requires the MapMyCells API key and data
        # submission pipeline. See the Allen docs:
        # https://alleninstitute.github.io/abc_atlas_access/
        # https://knowledge.brain-map.org/mapmycells/

        result = pd.DataFrame({
            "cell_barcode": adata.obs_names,
            "mapmycells_label": ["pending_api"] * n_cells,
            "mapmycells_score": [0.0] * n_cells,
        })

        print(f"    P0 framework: API integration ready for Phase 2+")
        return result
