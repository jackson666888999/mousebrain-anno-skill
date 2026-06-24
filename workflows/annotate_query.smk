# Snakemake workflow: Annotate query data
# Run: snakemake -s workflows/annotate_query.smk --config input=query.h5ad ref=allen-consensus-wmb level=subclass

rule all:
    input:
        "results/query.annotated.h5ad",
        "results/query.annotations.tsv",
        "results/query.method_report.html",

rule annotate:
    input:
        "data/references/{ref}_vlatest_metadata.json",
    output:
        "results/query.annotated.h5ad",
    params:
        ref=config.get("ref", "allen-consensus-wmb"),
        level=config.get("level", "subclass"),
        methods=config.get("methods", "mapmycells,scanvi,celltypist"),
        region=config.get("region", ""),
    shell:
        "mbanno annotate -i {config[input]} -o results/ "
        "--reference {params.ref} --level {params.level} "
        "--methods {params.methods} "
        "{'--region ' + params.region if params.region else ''}"
