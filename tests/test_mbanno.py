"""
mbanno — Unit tests

Run: pytest tests/ -v
"""

import sys
from pathlib import Path

# Ensure package is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np
import pandas as pd
import pytest


class TestImports:
    """Test that all core modules can be imported."""

    def test_import_mbanno(self):
        from mbanno import __version__, ConsensusAnnotator, QCProcessor
        assert __version__ == "0.1.0"
        assert ConsensusAnnotator is not None
        assert QCProcessor is not None

    def test_import_io(self):
        from mbanno import io
        assert hasattr(io, "read_adata")
        assert hasattr(io, "write_adata")

    def test_import_consensus(self):
        from mbanno.consensus import ConsensusAnnotator
        annotator = ConsensusAnnotator()
        assert annotator.reference == "allen-consensus-wmb"
        assert annotator.taxonomy_level == "subclass"

    def test_import_confidence(self):
        from mbanno.confidence import ConfidenceEvaluator
        evaluator = ConfidenceEvaluator()
        assert evaluator.min_confidence == 0.8

    def test_import_qc(self):
        from mbanno.qc import QCProcessor
        qc = QCProcessor(min_genes=200, min_counts=500, max_pct_mt=20.0)
        assert qc.min_genes == 200
        assert qc.min_counts == 500

    def test_import_marker_rules(self):
        from mbanno.marker_rules import MarkerRuleEngine
        engine = MarkerRuleEngine()
        rules = engine._get_default_rules()
        assert "neuron" in rules
        assert "glia" in rules

    def test_import_report(self):
        from mbanno.report import ReportGenerator
        assert ReportGenerator is not None

    def test_import_benchmark(self):
        from mbanno.benchmark import BenchmarkRunner
        assert BenchmarkRunner is not None

    def test_import_method_runners(self):
        from mbanno.mapmycells import MapMyCellsRunner
        from mbanno.scanvi import ScanviAnnotator
        from mbanno.celltypist_runner import CelltypistAnnotator
        from mbanno.singler_runner import SinglerAnnotator
        assert MapMyCellsRunner is not None
        assert ScanviAnnotator is not None
        assert CelltypistAnnotator is not None
        assert SinglerAnnotator is not None


class TestConfidence:
    """Test confidence evaluation logic."""

    def test_confidence_computation(self):
        from mbanno.confidence import ConfidenceEvaluator

        evaluator = ConfidenceEvaluator(min_confidence=0.5)
        n = 100

        annotations = pd.DataFrame({
            "cell_barcode": [f"cell_{i}" for i in range(n)],
            "mapmycells_label": ["type_A"] * 50 + ["type_B"] * 50,
            "scanvi_label": (["type_A"] * 45 + ["type_B"] * 5 +
                             ["type_A"] * 5 + ["type_B"] * 45),
            "mapmycells_score": np.random.uniform(0.5, 1.0, n),
            "scanvi_score": np.random.uniform(0.5, 1.0, n),
        })

        result = evaluator.evaluate(None, annotations)

        assert "confidence" in result.columns
        assert "confidence_level" in result.columns
        assert "agreement_rate" in result.columns
        assert len(result) == n

        # High agreement cells should have higher confidence
        high_agreement = result[result["agreement_rate"] >= 0.8]
        low_agreement = result[result["agreement_rate"] < 0.8]
        if len(high_agreement) > 0 and len(low_agreement) > 0:
            assert high_agreement["confidence"].mean() > low_agreement["confidence"].mean()

    def test_confidence_levels(self):
        from mbanno.confidence import ConfidenceEvaluator

        evaluator = ConfidenceEvaluator(min_confidence=0.5)
        n = 10

        annotations = pd.DataFrame({
            "cell_barcode": [f"cell_{i}" for i in range(n)],
            "mapmycells_label": ["type_A"] * n,
            "scanvi_label": ["type_A"] * n,
            "mapmycells_score": [0.9] * n,
            "scanvi_score": [0.9] * n,
        })

        result = evaluator.evaluate(None, annotations)

        valid_levels = {"high_confidence", "medium_confidence", "low_confidence", "ambiguous", "novel_or_state"}
        for level in result["confidence_level"]:
            assert level in valid_levels


class TestMarkerRules:
    """Test marker rule engine."""

    def test_default_rules_load(self):
        from mbanno.marker_rules import MarkerRuleEngine
        engine = MarkerRuleEngine()
        rules = engine._get_default_rules()

        assert "glutamatergic" in rules.get("neuron", {})
        assert "astrocyte" in rules.get("glia", {})
        assert "microglia" in rules.get("glia", {})
        assert "endothelial" in rules.get("vascular", {})

    def test_glutamatergic_markers(self):
        from mbanno.marker_rules import MarkerRuleEngine
        engine = MarkerRuleEngine()
        rules = engine._get_default_rules()

        glut = rules["neuron"]["glutamatergic"]
        assert "Slc17a7" in glut  # VGLUT1
        assert "Camk2a" in glut

    def test_GABAergic_markers(self):
        from mbanno.marker_rules import MarkerRuleEngine
        engine = MarkerRuleEngine()
        rules = engine._get_default_rules()

        gaba = rules["neuron"]["GABAergic"]
        assert "Gad1" in gaba
        assert "Gad2" in gaba
        assert "Slc32a1" in gaba  # VGAT


class TestBenchmark:
    """Test benchmark framework."""

    def test_metric_computation(self):
        from mbanno.benchmark import BenchmarkRunner

        runner = BenchmarkRunner.__new__(BenchmarkRunner)
        y_true = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        y_pred = np.array([0, 0, 1, 1, 1, 2, 2, 2, 2])

        metrics = runner._compute_metrics(y_true, y_pred)
        assert "accuracy" in metrics
        assert "macro_f1" in metrics
        assert "weighted_f1" in metrics
        assert "balanced_accuracy" in metrics
        assert "ari" in metrics
        assert 0 <= metrics["accuracy"] <= 1


class TestQC:
    """Test QC processor."""

    def test_qc_init_params(self):
        from mbanno.qc import QCProcessor

        qc = QCProcessor(min_genes=100, min_counts=300, max_pct_mt=15.0)
        assert qc.min_genes == 100
        assert qc.min_counts == 300
        assert qc.max_pct_mt == 15.0
        assert qc.doublet_threshold == 0.25

    def test_qc_summary_empty(self):
        from mbanno.qc import QCProcessor
        import anndata

        qc = QCProcessor()
        adata = anndata.AnnData(X=np.random.poisson(5, (100, 50)))
        adata.obs["n_genes_by_counts"] = np.random.randint(50, 5000, 100)
        adata.obs["total_counts"] = np.random.randint(500, 50000, 100)
        adata.obs["pct_counts_mt"] = np.random.uniform(0, 15, 100)
        adata.obs["predicted_doublet"] = np.random.choice([True, False], 100, p=[0.05, 0.95])

        summary = qc.get_summary(adata)
        assert "n_cells" in summary
        assert "median_genes_per_cell" in summary
        assert "doublet_rate" in summary
