# Snakemake workflow: Download reference datasets
# Run: snakemake -s workflows/download_references.smk -j4

DATASETS = [
    "allen-consensus-wmb",
    "wmb-2023",
]

rule all:
    input:
        expand("data/references/{dataset}_vlatest_metadata.json", dataset=DATASETS),
        "data/references/markers/marker_genes.json",

rule download_reference:
    output:
        "data/references/{dataset}_vlatest_metadata.json",
    shell:
        "mbanno download-reference --ref {wildcards.dataset}"

rule download_markers:
    output:
        "data/references/markers/marker_genes.json",
    shell:
        "mbanno download-reference --markers-only"
