"""
Benchmark — Reproducible evaluation framework for cell type annotation

Evaluates annotation methods on held-out Allen/BICCN reference data
with standard metrics: F1, accuracy, hierarchical F1, ARI, etc.

Benchmark claims are only made when supported by reproducible results.
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import (
    accuracy_score,
    adjusted_rand_score,
    balanced_accuracy_score,
    f1_score,
)


class BenchmarkRunner:
    """Run reproducible benchmarks on annotation methods.

    Parameters
    ----------
    config_path : str
        Path to YAML benchmark configuration.
    output_dir : str
        Directory for benchmark results.
    """

    def __init__(self, config_path: str, output_dir: str = "benchmark_results"):
        self.config_path = Path(config_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        with open(self.config_path) as f:
            self.config = yaml.safe_load(f)

        self.results = []

    def run(self, methods: List[str]):
        """Run benchmarks for specified methods."""
        print(f"\n{'=' * 60}")
        print(f"  mbanno Benchmark")
        print(f"{'=' * 60}")
        print(f"  Config:    {self.config_path}")
        print(f"  Methods:   {', '.join(methods)}")
        print(f"  Output:    {self.output_dir}")
        print(f"{'=' * 60}\n")

        for task in self.config.get("tasks", []):
            print(f"  Task: {task['name']}")
            result = self._run_task(task, methods)
            self.results.append(result)

        self._save_results()
        self._print_summary()

    def _run_task(self, task: Dict, methods: List[str]) -> Dict:
        """Run a single benchmark task."""
        task_name = task["name"]
        dataset = task.get("dataset", "")
        split = task.get("split", "random")
        level = task.get("level", "subclass")

        metrics = {}
        for method in methods:
            # Placeholder: In Phase 2+, load actual predictions
            # and compare to gold-standard labels
            metrics[method] = {
                "accuracy": 0.0,
                "macro_f1": 0.0,
                "weighted_f1": 0.0,
                "balanced_accuracy": 0.0,
                "ari": 0.0,
            }

        return {
            "task": task_name,
            "dataset": dataset,
            "split": split,
            "level": level,
            "methods": metrics,
            "timestamp": datetime.now().isoformat(),
        }

    def _compute_metrics(self, y_true: np.ndarray, y_pred: np.ndarray) -> Dict:
        """Compute standard annotation metrics."""
        labels = np.unique(np.concatenate([y_true, y_pred]))
        return {
            "accuracy": float(accuracy_score(y_true, y_pred)),
            "macro_f1": float(f1_score(y_true, y_pred, average="macro", labels=labels, zero_division=0)),
            "weighted_f1": float(f1_score(y_true, y_pred, average="weighted", labels=labels, zero_division=0)),
            "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
            "ari": float(adjusted_rand_score(y_true, y_pred)),
        }

    def _save_results(self):
        """Save benchmark results to JSON and CSV."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # JSON
        json_path = self.output_dir / f"benchmark_{timestamp}.json"
        with open(json_path, "w") as f:
            json.dump(self.results, f, indent=2, default=str)
        print(f"\n  Results saved: {json_path}")

        # Summary CSV
        rows = []
        for task in self.results:
            for method, scores in task["methods"].items():
                rows.append({
                    "task": task["task"],
                    "level": task["level"],
                    "method": method,
                    **scores,
                })
        if rows:
            df = pd.DataFrame(rows)
            csv_path = self.output_dir / f"benchmark_{timestamp}.csv"
            df.to_csv(csv_path, index=False)
            print(f"  Summary saved: {csv_path}")

    def _print_summary(self):
        """Print a summary table of benchmark results."""
        print(f"\n{'=' * 60}")
        print(f"  Benchmark Summary")
        print(f"{'=' * 60}")
        print(f"  Tasks run: {len(self.results)}")
        print(f"  Results in: {self.output_dir}/")
        print(f"  Data sources: {self.config.get('data_sources', 'data_sources.yaml')}")
        print(f"{'=' * 60}")
        print(f"\n  ⚠ Phase 1: Benchmark framework ready.")
        print(f"  Full results after Phase 2+ integration with")
        print(f"  held-out Allen/BICCN reference data.")
        print()
