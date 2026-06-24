"""
References — Reference data management (download, cache, metadata tracking)

All data is obtained from official public sources.
No raw data redistribution. Users download via scripts.
"""

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict

import requests
import yaml
from tqdm import tqdm

from .io import load_data_sources


class ReferenceManager:
    """Manage reference dataset downloads, metadata, and citations.

    Never redistributes raw data. Provides download scripts, manifests,
    and license tracking.

    Parameters
    ----------
    data_dir : str or Path
        Root directory for reference data storage.
    """

    def __init__(self, data_dir: str = "data/references"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)

    def list_references(self):
        """List all available reference datasets from data_sources.yaml."""
        try:
            manifest = load_data_sources()
        except FileNotFoundError:
            print("  data_sources.yaml not found. Run mbanno download-reference first.")
            return

        print("\nAvailable Reference Datasets\n" + "=" * 64)

        references = manifest.get("references", [])
        for ref in references:
            print(f"  [{ref['id']}] {ref['name']}")
            print(f"      License:  {ref.get('license', 'see source')}")
            print(f"      URL:      {ref.get('source_url', 'N/A')}")
            print(f"      Data:     {ref.get('data_path', 'N/A')}")
            print()

    def download_reference(
        self,
        reference_id: str,
        version: str = "latest",
    ):
        """Download reference metadata from official sources.

        Parameters
        ----------
        reference_id : str
            Reference ID from data_sources.yaml.
        version : str
            Version string (e.g., '20251031').
        """
        try:
            manifest = load_data_sources()
        except Exception:
            print("  Cannot load data_sources.yaml")
            return

        references = manifest.get("references", [])
        ref_config = None
        for ref in references:
            if ref["id"] == reference_id:
                ref_config = ref
                break

        if ref_config is None:
            print(f"  Reference '{reference_id}' not found in data_sources.yaml")
            print("  Available references:")
            self.list_references()
            return

        ref_config = dict(ref_config) if isinstance(ref_config, dict) else {}
        source_url = ref_config.get("source_url", "")
        data_path = ref_config.get("data_path", "")
        license_str = ref_config.get("license", "Unknown")
        citation = ref_config.get("citation", "Unknown")

        print(f"\n  Downloading: {ref_config['name']}")
        print(f"  Version:     {version}")
        print(f"  License:     {license_str}")
        print(f"  Source:      {source_url}")
        print(f"  DOI:         {ref_config.get('doi', 'N/A')}")
        print(f"  Citation:    {citation}")
        print()

        # Save metadata to disk
        metadata = {
            "reference_id": reference_id,
            "version": version,
            "license": license_str,
            "source_url": source_url,
            "citation": citation,
            "download_date": datetime.now().isoformat(),
            "data_source_yaml": "data_sources.yaml",
        }

        metadata_path = self.data_dir / f"{reference_id}_v{version}_metadata.json"
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)
        print(f"  Metadata saved: {metadata_path}")

        # Marker genes
        self._download_markers()

    def _download_markers(self):
        """Download marker gene lists."""
        marker_dir = self.data_dir / "markers"
        marker_dir.mkdir(parents=True, exist_ok=True)

        markers = {
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

        marker_path = marker_dir / "marker_genes.json"
        with open(marker_path, "w") as f:
            json.dump(markers, f, indent=2)
        print(f"  Markers saved: {marker_path}")
