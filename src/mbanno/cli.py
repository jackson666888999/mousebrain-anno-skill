"""
CLI — mbanno command-line interface
"""

import click
import sys
from pathlib import Path
from typing import Optional


@click.group()
@click.version_option(version="0.1.0", prog_name="mbanno")
def main():
    """mbanno — Mouse Brain Consensus Annotator

    A reproducible consensus annotation framework for mouse whole-brain
    single-cell and spatial transcriptomics using Allen/BICCN reference taxonomies.
    """
    pass


@main.command()
@click.option("--ref", "-r", "reference", default="allen-consensus-wmb",
              help="Reference dataset ID")
@click.option("--version", "-v", default="20251031",
              help="Reference version")
@click.option("--output-dir", "-o", default="data/references",
              help="Output directory")
@click.option("--list", "list_refs", is_flag=True, help="List available references")
def download_reference(reference: str, version: str, output_dir: str, list_refs: bool):
    """Download reference data from official sources.

    Downloads metadata, taxonomy files, and marker gene lists from
    Allen Institute / BICCN public repositories.

    Examples:

    \b
        mbanno download-reference --ref allen-consensus-wmb --version 20251031
        mbanno download-reference --list
    """
    from .references import ReferenceManager

    if list_refs:
        rm = ReferenceManager()
        rm.list_references()
        return

    rm = ReferenceManager(data_dir=output_dir)
    rm.download_reference(reference_id=reference, version=version)


@main.command()
@click.option("--input", "-i", "input_path", required=True,
              help="Input h5ad file path")
@click.option("--output", "-o", "output_dir", default="results",
              help="Output directory")
@click.option("--reference", "-r", default="allen-consensus-wmb",
              help="Reference dataset ID")
@click.option("--methods", "-m", multiple=True,
              default=["mapmycells", "scanvi", "celltypist"],
              help="Annotation methods to use")
@click.option("--level", "-l", default="subclass",
              type=click.Choice(["class", "subclass", "supertype", "cluster"]),
              help="Taxonomy level for annotation")
@click.option("--region", default=None,
              help="Brain region constraint (e.g., cortex, hippocampus)")
@click.option("--min-confidence", default=0.5, type=float,
              help="Minimum confidence threshold")
@click.option("--no-report", is_flag=True, help="Skip HTML report generation")
def annotate(
    input_path: str,
    output_dir: str,
    reference: str,
    methods: tuple,
    level: str,
    region: Optional[str],
    min_confidence: float,
    no_report: bool,
):
    """Annotate mouse brain single-cell data.

    Runs multi-tool consensus annotation on an h5ad file using specified
    methods and reference taxonomies.

    Examples:

    \b
        mbanno annotate -i query.h5ad -o results/
        mbanno annotate -i query.h5ad -m mapmycells scanvi celltypist marker
        mbanno annotate -i query.h5ad --region cortex --level subclass
    """
    import scanpy as sc
    from .consensus import ConsensusAnnotator

    methods_list = list(methods)

    print(f"\n{'=' * 60}")
    print(f"mbanno — Annotation Pipeline")
    print(f"{'=' * 60}")
    print(f"  Input:       {input_path}")
    print(f"  Reference:   {reference}")
    print(f"  Methods:     {', '.join(methods_list)}")
    print(f"  Level:       {level}")
    print(f"  Region:      {region or 'none'}")
    print(f"{'=' * 60}\n")

    adata = sc.read_h5ad(input_path)
    print(f"  Cells: {adata.n_obs:,}  Genes: {adata.n_vars:,}")

    annotator = ConsensusAnnotator(
        reference=reference,
        taxonomy_level=level,
        data_dir="data/references",
    )

    results = annotator.annotate(
        adata=adata,
        methods=methods_list,
        region_constraint=region,
        min_confidence=min_confidence,
    )

    annotator.save_results(results, output_dir=Path(output_dir))

    if not no_report:
        annotator.generate_report(results, output_dir=Path(output_dir))

    print(f"\n  Done. Results saved to: {output_dir}/\n")


@main.command()
@click.option("--config", "-c", "config_path", required=True,
              help="Benchmark configuration YAML file")
@click.option("--output", "-o", "output_dir", default="benchmark_results",
              help="Output directory")
@click.option("--methods", "-m", multiple=True,
              default=["mapmycells", "scanvi", "celltypist", "singler"],
              help="Methods to benchmark")
def benchmark(config_path: str, output_dir: str, methods: tuple):
    """Run benchmark evaluation.

    Evaluates annotation methods on held-out public datasets using
    standard metrics (F1, accuracy, hierarchical F1, etc.).

    Examples:

    \b
        mbanno benchmark -c benchmark.yaml -o results/
    """
    from .benchmark import BenchmarkRunner

    runner = BenchmarkRunner(config_path=config_path, output_dir=output_dir)
    runner.run(methods=list(methods))


@main.command()
@click.option("--input", "-i", "input_path", required=True,
              help="Input annotated h5ad file")
@click.option("--output", "-o", "output_dir", default="results",
              help="Output directory")
@click.option("--format", "-f", default="html",
              type=click.Choice(["html", "markdown", "json"]),
              help="Report format")
def report(input_path: str, output_dir: str, format: str):
    """Generate annotation report from results.

    Creates a comprehensive HTML/Markdown/JSON report with UMAP
    visualizations, confidence distributions, and method agreement.

    Examples:

    \b
        mbanno report -i results/query.annotated.h5ad
        mbanno report -i results/query.annotated.h5ad -f markdown
    """
    import scanpy as sc
    from .report import ReportGenerator

    adata = sc.read_h5ad(input_path)
    reporter = ReportGenerator(output_dir=output_dir)
    reporter.generate(adata, format=format)
    print(f"\n  Report saved to: {output_dir}/\n")


if __name__ == "__main__":
    main()
