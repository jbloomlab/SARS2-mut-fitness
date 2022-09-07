"""Top-level ``snakemake`` file that runs pipeline."""


import glob


configfile: "config.yaml"


rule all:
    """Target rule with desired output files."""
    input:
        "_temp.txt"


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


checkpoint samples_by_clade:
    """Get samples in mutation-annotated tree by nextstrain clade."""
    input:
        csv=rules.mat_samples.output.csv,
    output:
        subdir=directory("results/mat/samples_by_clade"),
    script:
        "scripts/samples_by_clade.py"


def clades(wc):
    """Return list of all clades."""
    subdir = checkpoints.samples_by_clade.get(**wc).output.subdir
    return [
        os.path.splitext(os.path.basename(f))[0] for f in glob.glob(f"{subdir}/*.txt")
    ]


rule mat_clade_subset:
    """Get mutation-annotated tree for just a clade."""
    input:
        mat=rules.get_mat_tree.output.mat,
        samples=os.path.join(rules.samples_by_clade.output.subdir, "{clade}.txt"),
    output:
        mat="results/mat/mats_by_clade/{clade}_mat_tree.pb",
    shell:
        "matUtils extract -i {input.mat} -s {input.samples} -o {output.mat}"


rule _temp:
    input:
        lambda wc: [
            f"results/mat/mats_by_clade/{clade}_mat_tree.pb" for clade in clades(wc)
        ],
    output:
        "_temp.txt"
    shell:
        "echo not_implemented"
