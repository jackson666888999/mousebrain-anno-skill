"""
Confidence — Multi-dimensional confidence evaluation for cell type annotations

Evaluates: agreement rate, prediction margin, marker support,
region consistency, and overall confidence level.
"""

from typing import Dict, Optional

import anndata
import numpy as np
import pandas as pd
import scipy.stats


class ConfidenceEvaluator:
    """Evaluate annotation confidence using multiple evidence sources.

    Produces per-cell confidence scores and categorical confidence levels:
    high_confidence, medium_confidence, low_confidence, ambiguous, novel_or_state.

    Parameters
    ----------
    min_confidence : float
        Threshold for high-confidence calls (default 0.8).
    """

    def __init__(self, min_confidence: float = 0.8):
        self.min_confidence = min_confidence

    def evaluate(
        self,
        adata: anndata.AnnData,
        annotations: pd.DataFrame,
        min_confidence: Optional[float] = None,
    ) -> pd.DataFrame:
        """Run full confidence evaluation.

        Parameters
        ----------
        adata : AnnData
            Expression data.
        annotations : pd.DataFrame
            Annotation results with per-method labels.
        min_confidence : float, optional
            Override default confidence threshold.

        Returns
        -------
        pd.DataFrame
            Per-cell confidence metrics.
        """
        min_conf = min_confidence or self.min_confidence

        results = []

        n_cells = len(annotations)

        # 1. Agreement rate (fraction of methods agreeing)
        label_cols = [c for c in annotations.columns if c.endswith("_label")]
        if len(label_cols) >= 2:
            label_matrix = annotations[label_cols].values
            agreement = np.array([
                (label_matrix[i] == label_matrix[i][0]).mean()
                if len(label_matrix[i]) > 0 else 0.0
                for i in range(n_cells)
            ])
        else:
            agreement = np.ones(n_cells)

        # 2. Prediction margin (if score columns exist)
        score_cols = [c for c in annotations.columns if c.endswith("_score")]
        if score_cols:
            score_matrix = annotations[score_cols].values
            sorted_scores = np.sort(score_matrix, axis=1)[:, ::-1]
            if score_matrix.shape[1] >= 2:
                margins = sorted_scores[:, 0] - sorted_scores[:, 1]
            else:
                margins = sorted_scores[:, 0]
            margins = (margins - margins.min()) / (margins.max() - margins.min() + 1e-8)
        else:
            margins = np.ones(n_cells) * 0.5

        # 3. Overall confidence (weighted combination of agreement + margin)
        confidence = 0.6 * agreement + 0.4 * margins

        # 4. Categorical confidence levels
        levels = np.full(n_cells, "medium_confidence", dtype=object)
        levels[confidence >= 0.8] = "high_confidence"
        levels[confidence < min_conf] = "low_confidence"
        # If only one method and low margin, ambiguous
        if len(label_cols) == 1:
            levels[(confidence < 0.5)] = "ambiguous"

        # 5. Estimated entropy
        entropy = np.ones(n_cells) * 0.5
        if len(label_cols) >= 2:
            for i in range(n_cells):
                vals, counts = np.unique(label_matrix[i], return_counts=True)
                probs = counts / counts.sum()
                entropy[i] = float(scipy.stats.entropy(probs))

        result = pd.DataFrame({
            "cell_barcode": annotations.get("cell_barcode", np.arange(n_cells)),
            "agreement_rate": agreement,
            "prediction_margin": margins,
            "confidence": confidence,
            "entropy": entropy,
            "confidence_level": levels,
            "marker_support": np.zeros(n_cells),  # placeholder for Phase 2 integration
            "region_consistency": np.ones(n_cells),  # placeholder
        })

        return result
