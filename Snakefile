"""Top-level ``snakemake`` file that runs pipeline."""


import glob

import pandas as pd


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


checkpoint mat_samples:
    """Get all samples in mutation-annotated tree with their dates and clades."""
    input:
        mat=rules.get_mat_tree.output.mat,
    output:
        csv="results/mat/samples.csv",
        clade_counts="results/mat/sample_clade_counts.csv",
    params:
        min_clade_samples=config["min_clade_samples"],
    script:
        "scripts/mat_samples.py"


def clades_w_adequate_counts(wc):
    """Return list of all clades with adequate sample counts."""
    return (
        pd.read_csv(checkpoints.mat_samples.get(**wc).output.clade_counts)
        .query("adequate_sample_counts")
        ["nextstrain_clade"]
        .tolist()
    )


rule samples_by_clade:
    """Get samples in mutation-annotated tree by nextstrain clade."""
    input:
        csv=rules.mat_samples.output.csv,
    output:
        txt="results/mat_by_clade_subset/{clade}.txt",
    run:
        (
            pd.read_csv(input.csv)
            .query("nextstrain_clade == @wildcards.clade")
            ["sample"]
            .to_csv(output.txt, index=False, header=False)
        )


rule mat_clade_subset:
    """Get mutation-annotated tree for just a clade."""
    input:
        mat=rules.get_mat_tree.output.mat,
        samples=rules.samples_by_clade.output.txt,
    output:
        mat="results/mat_by_clade_subset/{clade}_mat_tree.pb",
    shell:
        "matUtils extract -i {input.mat} -s {input.samples} -o {output.mat}"


rule translate_mat:
    """Translate mutations on mutation-annotated tree for clade."""
    input:
        mat=rules.mat_clade_subset.output.mat,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
        ref_gtf=rules.get_ref_gtf.output.ref_gtf,
    output:
        tsv="results/mat_by_clade_subset/{clade}_translated_mutations.tsv",
    shell:
        """
        matUtils summary \
            -i {input.mat} \
            -g {input.ref_gtf} \
            -f {input.ref_fasta} \
            -t {output.tsv}
        """


rule clade_founder_json:
    """Get JSON with nexstrain clade founders (indels not included)."""
    params:
        url=config["clade_founder_json"],
    output:
        json="results/clade_founders_no_indels/clade_founders.json",
    shell:
        "curl {params.url} > {output.json}"


rule clade_founder_fasta:
    """Get FASTA for a nextstrain clade founder (indels not included)."""
    input:
        json=rules.clade_founder_json.output.json,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
    output:
        fasta="results/clade_founders_no_indels/{clade}.fa",
    script:
        "scripts/clade_founder_fasta.py"


rule count_mutations:
    """Count mutations, excluding branches with too many mutations or reversions."""
    input:
        tsv=rules.translate_mat.output.tsv,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
        clade_founder_fasta=rules.clade_founder_fasta.output.fasta,
    output:
        csv="results/mutation_counts/counts_by_clade/{clade}.csv",
    params:
        max_nt_mutations=config["max_nt_mutations"],
        max_reversions_to_ref=config["max_reversions_to_ref"],
        max_reversions_to_clade_founder=config["max_reversions_to_clade_founder"],
    log:
        notebook="results/mutation_counts/counts_by_clade/{clade}_count_mutations.ipynb",
    notebook:
        "notebooks/count_mutations.py.ipynb"


rule synonymous_mut_rates:
    """Compute overall rates of synonymous mutations."""
    input:
        counts=lambda wc: [
            f"results/mutation_counts/counts_by_clade/{clade}.csv"
            for clade in clades_w_adequate_counts(wc)
        ],
    output:
        "_temp.txt",
        csv="results/synonymous_mut_rates/rates.csv",
    log:
        notebook="results/synonymous_mut_rates/synonymous_mut_rates.ipynb",
    notebook:
        "notebooks/synonymous_mut_rates.py.ipynb"
