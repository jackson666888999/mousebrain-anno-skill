# Snakemake workflow: Run public benchmarks
# Run: snakemake -s workflows/benchmark_public.smk -j4

rule all:
    input:
        "benchmark_results/benchmark_summary.csv",

rule benchmark:
    input:
        "benchmark.yaml",
    output:
        "benchmark_results/benchmark_summary.csv",
    shell:
        "mbanno benchmark -c {input} -o benchmark_results/"
