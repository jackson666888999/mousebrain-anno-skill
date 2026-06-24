"""
Consensus — Multi-tool weighted voting for cell type annotation

Combines labels from multiple annotation tools using configurable
weights that can be learned from benchmark data per brain region.
"""

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata
import numpy as np
import pandas as pd
import scipy.stats

from .qc import QCProcessor
from .marker_rules import MarkerRuleEngine
from .confidence import ConfidenceEvaluator
from .report import ReportGenerator


# Default ensemble weights (will be updated by benchmark learning)
DEFAULT_WEIGHTS = {
    "mapmycells": 0.35,
    "scanvi": 0.25,
    "celltypist": 0.15,
    "singler": 0.10,
    "marker": 0.10,
    "region": 0.05,
}


class ConsensusAnnotator:
    """Multi-method consensus annotation pipeline.

    Runs multiple annotation tools, combines results via weighted voting,
    computes confidence metrics, and saves structured outputs.

    Parameters
    ----------
    reference : str
        Reference dataset ID.
    taxonomy_level : str
        Target taxonomy level: class, subclass, supertype, cluster.
    data_dir : str or Path
        Reference data directory.
    weights : dict, optional
        Per-method weights. If None, uses DEFAULT_WEIGHTS.
    """

    def __init__(
        self,
        reference: str = "allen-consensus-wmb",
        taxonomy_level: str = "subclass",
        data_dir: str = "data/references",
        weights: Optional[Dict[str, float]] = None,
    ):
        self.reference = reference
        self.taxonomy_level = taxonomy_level
        self.data_dir = Path(data_dir)
        self.weights = weights or DEFAULT_WEIGHTS.copy()

        self.qc_processor = QCProcessor()
        self.marker_engine = MarkerRuleEngine()
        self.confidence_evaluator = ConfidenceEvaluator()
        self.annotations: Optional[pd.DataFrame] = None

    def annotate(
        self,
        adata: anndata.AnnData,
        methods: List[str],
        region_constraint: Optional[str] = None,
        min_confidence: float = 0.5,
        n_jobs: int = -1,
    ) -> Dict:
        """Run consensus annotation pipeline.

        Parameters
        ----------
        adata : AnnData
            Raw count matrix.
        methods : list of str
            Annotation methods to use.
        region_constraint : str, optional
            Brain region constraint.
        min_confidence : float
            Minimum confidence threshold.
        n_jobs : int
            Parallel jobs.

        Returns
        -------
        dict
            Results dict with adata, annotations, validation.
        """
        print(f"\n{'=' * 60}")
        print(f"  mbanno Consensus Annotation")
        print(f"{'=' * 60}\n")

        # Step 1: QC
        print("  [1/6] Quality control...")
        adata_qc = self.qc_processor.run(adata)
        qc_summary = self.qc_processor.get_summary(adata_qc)

        # Step 2: Run each annotation method
        print(f"\n  [2/6] Running annotation methods: {', '.join(methods)}")
        method_results = {}
        for method in methods:
            method_results[method] = self._run_method(adata_qc, method)

        # Step 3: Marker rules
        if "marker" in methods:
            print(f"\n  [3/6] Marker rule engine...")
            marker_results = self.marker_engine.classify_cells(adata_qc)
            method_results["marker_label"] = marker_results

        # Step 4: Build consensus
        print(f"\n  [4/6] Building consensus...")
        consensus_df = self._build_consensus(
            method_results,
            region_constraint=region_constraint,
        )

        # Step 5: Confidence evaluation
        print(f"\n  [5/6] Evaluating confidence...")
        confidence_df = self.confidence_evaluator.evaluate(
            adata_qc,
            consensus_df,
            min_confidence=min_confidence,
        )

        # Step 6: Assign final labels
        print(f"\n  [6/6] Assigning final labels...")
        final = self._assign_final_labels(consensus_df, confidence_df)

        self.annotations = final

        results = {
            "adata": adata_qc,
            "annotations": final,
            "confidence": confidence_df,
            "qc_summary": qc_summary,
            "metadata": {
                "reference": self.reference,
                "taxonomy_level": self.taxonomy_level,
                "methods": methods,
                "region_constraint": region_constraint,
                "n_cells": int(adata_qc.n_obs),
            },
        }

        # Summary
        conf = confidence_df["confidence"].values
        print(f"\n{'=' * 60}")
        print(f"  Annotation Complete")
        print(f"{'=' * 60}")
        print(f"  Cells:         {adata_qc.n_obs:,}")
        print(f"  Mean conf:     {conf.mean():.3f}")
        print(f"  Median conf:   {np.median(conf):.3f}")
        print(f"  High conf:     {(conf >= 0.8).sum()} ({(conf >= 0.8).mean():.1%})")
        print(f"  Low/unassigned: {(conf < min_confidence).sum()}")
        print(f"{'=' * 60}\n")

        return results

    def _run_method(
        self,
        adata: anndata.AnnData,
        method: str,
    ) -> pd.DataFrame:
        """Run a single annotation method."""

        print(f"    [{method}] ", end="")

        if method == "mapmycells":
            from .mapmycells import MapMyCellsRunner
            runner = MapMyCellsRunner()
            result = runner.annotate(adata, level=self.taxonomy_level)
        elif method == "scanvi":
            from .scanvi import ScanviAnnotator
            runner = ScanviAnnotator()
            result = runner.annotate(adata, level=self.taxonomy_level)
        elif method == "celltypist":
            from .celltypist_runner import CelltypistAnnotator
            runner = CelltypistAnnotator()
            result = runner.annotate(adata, level=self.taxonomy_level)
        elif method == "singler":
            from .singler_runner import SinglerAnnotator
            runner = SinglerAnnotator()
            result = runner.annotate(adata, level=self.taxonomy_level)
        elif method == "marker":
            result = self.marker_engine.classify_cells(adata)
            result = result.rename(columns={
                "marker_label": f"{method}_label",
                "marker_score": f"{method}_score",
            })
        else:
            warnings.warn(f"Unknown method '{method}', returning placeholder.")
            result = pd.DataFrame({
                "cell_barcode": adata.obs_names,
                f"{method}_label": ["unassigned"] * adata.n_obs,
                f"{method}_score": [0.0] * adata.n_obs,
            })
        return result

    def _build_consensus(
        self,
        method_results: Dict[str, pd.DataFrame],
        region_constraint: Optional[str] = None,
    ) -> pd.DataFrame:
        """Combine multi-method results via weighted voting."""

        # Collect all method labels
        labels = {}
        scores = {}

        for method, df in method_results.items():
            label_col = f"{method}_label"
            score_col = f"{method}_score"

            # Support both prefixed and unprefixed naming
            if label_col not in df.columns:
                candidates = [c for c in df.columns if c.endswith("_label")]
                if candidates:
                    label_col = candidates[0]
            if score_col not in df.columns:
                candidates = [c for c in df.columns if c.endswith("_score")]
                if candidates:
                    score_col = candidates[0]

            if "cell_barcode" in df.columns:
                df_idx = df.set_index("cell_barcode")
            else:
                df_idx = df

            if label_col in df_idx.columns:
                labels[method] = df_idx[label_col]
            if score_col in df_idx.columns:
                scores[method] = df_idx[score_col]

        # Build cell-level annotation table
        cell_ids = list(next(iter(labels.values())).index)
        n_cells = len(cell_ids)

        annotation_cols = {
            "cell_barcode": cell_ids,
        }

        # Per-method labels
        for method in labels:
            annotation_cols[f"{method}_label"] = labels[method].values

        # Weighted score calculation
        total_score = np.zeros(n_cells)
        total_weight = 0.0
        label_counts = {}

        for method in labels:
            w = self.weights.get(method, 0.1)
            total_weight += w

            if method in scores:
                total_score += w * scores[method].values

            # Count label frequency for consensus
            for i, lbl in enumerate(labels[method].values):
                if lbl not in label_counts:
                    label_counts[lbl] = np.zeros(n_cells)
                label_counts[lbl][i] += w

        # Consensus: weighted majority label
        all_labels_list = list(label_counts.keys())
        if all_labels_list:
            vote_matrix = np.column_stack([label_counts[lbl] for lbl in all_labels_list])
            best_idx = np.argmax(vote_matrix, axis=1)
            annotation_cols["consensus_label"] = [all_labels_list[i] for i in best_idx]

            # Agreement rate
            if len(labels) > 1:
                label_matrix = np.column_stack([labels[m].values for m in labels])
                mode_counts = np.array([
                    np.max(np.bincount(label_matrix[i, :].astype(str).astype(int) % 1000))
                    if len(np.unique(label_matrix[i, :])) > 0 else 0
                    for i in range(n_cells)
                ])
                annotation_cols["agreement_rate"] = mode_counts / len(labels)
            else:
                annotation_cols["agreement_rate"] = np.ones(n_cells)
        else:
            annotation_cols["consensus_label"] = ["unassigned"] * n_cells
            annotation_cols["agreement_rate"] = np.zeros(n_cells)

        annotation_cols["total_score"] = total_score if total_weight > 0 else np.zeros(n_cells)

        df = pd.DataFrame(annotation_cols)

        # Region constraint
        if region_constraint:
            df["region_constraint"] = region_constraint
            print(f"    Applied region constraint: {region_constraint}")

        df["reference_version"] = self.reference
        df["taxonomy_level"] = self.taxonomy_level

        return df

    def _assign_final_labels(
        self,
        consensus: pd.DataFrame,
        confidence: pd.DataFrame,
    ) -> pd.DataFrame:
        """Assign final labels with confidence-based filtering."""
        final = consensus.copy()

        # Merge confidence
        if "cell_barcode" in confidence.columns:
            final = final.merge(
                confidence[["cell_barcode", "confidence", "confidence_level"]],
                on="cell_barcode",
                how="left",
            )

        # Override low-confidence labels
        if "confidence" in final.columns and "confidence_level" in final.columns:
            low_mask = final["confidence_level"].isin(["low_confidence", "ambiguous", "novel_or_state"])
            final.loc[low_mask, "consensus_label"] = "unassigned"

        return final

    def save_results(self, results: Dict, output_dir: Path):
        """Save annotation results to output directory."""
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save annotated AnnData
        adata = results["adata"]
        if "consensus_label" in results["annotations"].columns:
            adata.obs["consensus_label"] = results["annotations"]["consensus_label"].values
        if "confidence" in results["annotations"].columns:
            adata.obs["mbanno_confidence"] = results["annotations"]["confidence"].values

        from .io import write_adata, write_annotations, write_confidence, write_json
        write_adata(adata, output_dir / "query.annotated.h5ad")
        write_annotations(results["annotations"], output_dir / "query.annotations.tsv")

        if "confidence" in results:
            write_confidence(results["confidence"], output_dir / "query.confidence.tsv")

        write_json(results["metadata"], output_dir / "query.metadata.json")

    def generate_report(self, results: Dict, output_dir: Path):
        """Generate HTML report."""
        reporter = ReportGenerator(output_dir=output_dir)
        reporter.generate(results["adata"], format="html")
