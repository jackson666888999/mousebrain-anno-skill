"""
Report — HTML/Markdown report generation for annotation results
"""

import datetime
from pathlib import Path
from typing import Union

import anndata
import numpy as np
import pandas as pd


class ReportGenerator:
    """Generate annotation reports in HTML or Markdown format.

    Parameters
    ----------
    output_dir : str or Path
        Report output directory.
    """

    def __init__(self, output_dir: Union[str, Path] = "results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate(
        self,
        adata: anndata.AnnData,
        format: str = "html",
    ):
        """Generate full annotation report.

        Parameters
        ----------
        adata : AnnData
            Annotated data object.
        format : str
            Report format: 'html', 'markdown', or 'json'.
        """
        print(f"  Generating report ({format})...")

        # Collect data
        summary = self._collect_summary(adata)

        if format == "html":
            self._generate_html(summary)
        elif format == "markdown":
            self._generate_markdown(summary)
        elif format == "json":
            self._generate_json(summary)

    def _collect_summary(self, adata: anndata.AnnData) -> dict:
        """Collect summary statistics from annotated data."""
        summary = {
            "timestamp": datetime.datetime.now().isoformat(),
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        }

        # Cell type distribution
        label_cols = [
            "consensus_label",
            "mapmycells_label",
            "scanvi_label",
            "celltypist_label",
            "singler_label",
            "marker_label",
        ]
        for col in label_cols:
            if col in adata.obs.columns:
                counts = adata.obs[col].value_counts().to_dict()
                summary[f"{col}_distribution"] = {
                    str(k): int(v) for k, v in counts.items()
                }

        # Confidence
        for col in ["mbanno_confidence", "confidence"]:
            if col in adata.obs.columns:
                conf = adata.obs[col].values
                summary["confidence"] = {
                    "mean": float(np.mean(conf)),
                    "median": float(np.median(conf)),
                    "std": float(np.std(conf)),
                    "high_conf_pct": float((conf >= 0.8).mean()),
                    "low_conf_pct": float((conf < 0.5).mean()),
                }
                break

        # QC metrics
        qc_cols = [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
        ]
        for col in qc_cols:
            if col in adata.obs.columns:
                vals = adata.obs[col].values
                summary[f"qc_{col}"] = {
                    "mean": float(np.mean(vals)),
                    "median": float(np.median(vals)),
                }

        return summary

    def _generate_html(self, summary: dict):
        """Generate an HTML report."""
        html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>mbanno Annotation Report</title>
<style>
body {{ font-family: Arial, sans-serif; max-width: 900px; margin: 40px auto; padding: 0 20px; }}
h1 {{ color: #333; border-bottom: 2px solid #44A2CD; padding-bottom: 10px; }}
h2 {{ color: #555; margin-top: 30px; }}
table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
th {{ background-color: #f4f4f4; }}
.stats {{ display: flex; gap: 20px; flex-wrap: wrap; }}
.card {{ background: #f9f9f9; border-left: 4px solid #44A2CD; padding: 15px 20px; border-radius: 4px; min-width: 150px; }}
.card h3 {{ margin: 0 0 5px 0; font-size: 12px; color: #888; text-transform: uppercase; }}
.card .value {{ font-size: 24px; font-weight: bold; color: #333; }}
.source {{ font-size: 11px; color: #999; margin-top: 40px; border-top: 1px solid #eee; padding-top: 20px; }}
</style>
</head>
<body>
<h1>mbanno — Mouse Brain Annotation Report</h1>
<p>Generated: {summary['timestamp']}</p>

<h2>Overview</h2>
<div class="stats">
<div class="card"><h3>Cells</h3><div class="value">{summary['n_cells']:,}</div></div>
<div class="card"><h3>Genes</h3><div class="value">{summary['n_genes']:,}</div></div>
"""

        # Confidence card
        if "confidence" in summary:
            c = summary["confidence"]
            html += f"""<div class="card"><h3>Mean Confidence</h3><div class="value">{c['mean']:.3f}</div></div>
<div class="card"><h3>High Conf</h3><div class="value">{c['high_conf_pct']:.1%}</div></div>
"""

        html += "</div>\n"

        # Cell type table
        html += "<h2>Cell Type Distribution</h2>\n"
        dist_key = None
        for key in summary:
            if key.endswith("_distribution"):
                dist_key = key
                break

        if dist_key:
            html += "<table><tr><th>Cell Type</th><th>Count</th><th>Percentage</th></tr>\n"
            dist = summary[dist_key]
            total = sum(dist.values())
            for ct, cnt in sorted(dist.items(), key=lambda x: -x[1]):
                pct = cnt / total * 100 if total > 0 else 0
                html += f"<tr><td>{ct}</td><td>{cnt:,}</td><td>{pct:.1f}%</td></tr>\n"
            html += "</table>\n"

        # QC table
        html += "<h2>QC Metrics</h2>\n"
        html += "<table><tr><th>Metric</th><th>Mean</th><th>Median</th></tr>\n"
        for key in sorted(summary):
            if key.startswith("qc_") and isinstance(summary[key], dict):
                metric_name = key.replace("qc_", "")
                html += f"<tr><td>{metric_name}</td><td>{summary[key]['mean']:.1f}</td><td>{summary[key]['median']:.1f}</td></tr>\n"
        html += "</table>\n"

        # Footer
        html += """<div class="source">
<p><strong>Data Sources:</strong> Allen/BICCN Whole Mouse Brain Atlas (Nature 2023, CC-BY-4.0).</p>
<p><strong>Citation:</strong> Yao, Z. et al. <em>Nature</em> 624, 317-332 (2023). doi:10.1038/s41586-023-06812-z</p>
<p><strong>Disclaimer:</strong> This report is for research use only. See <em>docs/data_licenses.md</em> for data use terms.</p>
</div>
</body>
</html>"""

        report_path = self.output_dir / "query.method_report.html"
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(html)
        print(f"  HTML report saved: {report_path}")

    def _generate_markdown(self, summary: dict):
        """Generate a Markdown report."""
        lines = [
            "# mbanno — Mouse Brain Annotation Report",
            f"Generated: {summary['timestamp']}",
            "",
            "## Overview",
            f"- Cells: {summary['n_cells']:,}",
            f"- Genes: {summary['n_genes']:,}",
        ]

        if "confidence" in summary:
            c = summary["confidence"]
            lines.extend([
                "",
                "## Confidence",
                f"- Mean: {c['mean']:.3f}",
                f"- Median: {c['median']:.3f}",
                f"- High conf (>0.8): {c['high_conf_pct']:.1%}",
                f"- Low conf (<0.5): {c['low_conf_pct']:.1%}",
            ])

        lines.extend([
            "",
            "## Data Sources",
            "Allen/BICCN Whole Mouse Brain Atlas (Nature 2023, CC-BY-4.0)",
            "Yao, Z. et al. *Nature* 624, 317-332 (2023)",
            "",
            "*Disclaimer: Research use only. See docs/data_licenses.md.*",
        ])

        report_path = self.output_dir / "query.method_report.md"
        with open(report_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))
        print(f"  Markdown report saved: {report_path}")

    def _generate_json(self, summary: dict):
        """Generate a JSON summary."""
        import json
        report_path = self.output_dir / "query.method_report.json"
        with open(report_path, "w") as f:
            json.dump(summary, f, indent=2, default=str)
        print(f"  JSON report saved: {report_path}")
