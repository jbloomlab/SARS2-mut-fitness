"""Top-level ``snakemake`` file that runs pipeline."""


configfile: "config.yaml"


rule all:
    """Target rule with desired output files."""
    input:
        "results/mat/samples_by_clade"


rule get_mat_tree:
    """Get the pre-built mutation-annotated tree."""
    params:
        url=config["mat_tree"],
    output:
        mat="results/mat/mat_tree.pb.gz"
    shell:
        "curl {params.url} > {output.mat}"


rule get_ref_fasta:
    """Get the reference FASTA."""
    params:
        url=config["ref_fasta"],
    output:
        ref_fasta="results/ref/ref.fa",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_fasta}"


rule get_ref_gtf:
    """Get the reference FASTA."""
    params:
        url=config["ref_gtf"],
    output:
        ref_gtf="results/ref/ref.gtf",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_gtf}"


rule mat_samples:
    """Get all samples in mutation-annotated tree with their dates and clades."""
    input:
        mat=rules.get_mat_tree.output.mat,
    output:
        csv="results/mat/samples.csv",
    script:
        "scripts/mat_samples.py"


rule mat_samples_by_clade:
    """Get samples in mutation-annotated tree by nextstrain clade."""
    input:
        csv=rules.mat_samples.output.csv,
    output:
        subdir=directory("results/mat/samples_by_clade"),
    script:
        "scripts/mat_samples_by_clade.py"
