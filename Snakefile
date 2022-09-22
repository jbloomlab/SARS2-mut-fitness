"""Top-level ``snakemake`` file that runs pipeline."""


import glob
import os

import Bio.SeqIO

import pandas as pd


configfile: "config.yaml"


rule all:
    """Target rule with desired output files."""
    input:
        "results/coding_nts/coding_nts.csv",


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


rule ref_coding_sites:
    """Get all sites in reference that are part of a coding sequence."""
    input:
        gtf=rules.get_ref_gtf.output.ref_gtf
    output:
        csv="results/ref/coding_sites.csv",
    script:
        "scripts/ref_coding_sites.py"


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


rule samples_by_clade_subset:
    """Get samples in mutation-annotated tree by nextstrain clade and subset."""
    input:
        csv=rules.mat_samples.output.csv,
    output:
        txt="results/mat_by_clade_subset/{clade}_{subset}.txt",
    params:
        match_regex=lambda wc: config["sample_subsets"][wc.subset]
    run:
        (
            pd.read_csv(input.csv)
            .query("nextstrain_clade == @wildcards.clade")
            .query(f"sample.str.match('{params.match_regex}')")
            ["sample"]
            .to_csv(output.txt, index=False, header=False)
        )


rule mat_clade_subset:
    """Get mutation-annotated tree for just a clade and subset."""
    input:
        mat=rules.get_mat_tree.output.mat,
        samples=rules.samples_by_clade_subset.output.txt,
    output:
        mat="results/mat_by_clade_subset/{clade}_{subset}.pb",
    shell:
        """
        if [ -s {input.samples} ]; then
            echo "Extracting samples from {input.samples}"
            matUtils extract -i {input.mat} -s {input.samples} -o {output.mat}
        else
            echo "No samples in {input.samples}"
            touch {output.mat}
        fi
        """


rule translate_mat:
    """Translate mutations on mutation-annotated tree for clade."""
    input:
        mat=rules.mat_clade_subset.output.mat,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
        ref_gtf=rules.get_ref_gtf.output.ref_gtf,
    output:
        tsv="results/mat_by_clade_subset/{clade}_{subset}_mutations.tsv",
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


rule clade_founder_fasta_and_muts:
    """Get FASTA and mutations for nextstrain clade founder (indels not included)."""
    input:
        json=rules.clade_founder_json.output.json,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
    output:
        fasta="results/clade_founders_no_indels/{clade}.fa",
        muts="results/clade_founders_no_indels/{clade}_ref_to_founder_muts.csv",
    script:
        "scripts/clade_founder_fasta.py"


rule count_mutations:
    """Count mutations, excluding branches with too many mutations or reversions."""
    input:
        tsv=rules.translate_mat.output.tsv,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
        clade_founder_fasta=rules.clade_founder_fasta_and_muts.output.fasta,
        ref_to_founder_muts=rules.clade_founder_fasta_and_muts.output.muts,
    output:
        csv="results/mutation_counts/{clade}_{subset}.csv",
    params:
        max_nt_mutations=config["max_nt_mutations"],
        max_reversions_to_ref=config["max_reversions_to_ref"],
        max_reversions_to_clade_founder=config["max_reversions_to_clade_founder"],
        exclude_ref_to_founder_muts=config["exclude_ref_to_founder_muts"],
        sites_to_exclude=config["sites_to_exclude"],
    log:
        notebook="results/mutation_counts/{clade}_{subset}_count_mutations.ipynb",
    notebook:
        "notebooks/count_mutations.py.ipynb"


rule aggregate_mutation_counts:
    """Aggregate the mutation counts for all clades and subsets."""
    input:
        counts=lambda wc: [
            f"results/mutation_counts/{clade}_{subset}.csv"
            for clade in clades_w_adequate_counts(wc)
            for subset in config["sample_subsets"]
        ],
    output:
        csv="results/mutation_counts/aggregated.csv",
    run:
        pd.concat(
            [
                pd.read_csv(f).assign(
                    clade=os.path.splitext(os.path.basename(f))[0].split("_")[0],
                    subset=os.path.splitext(os.path.basename(f))[0].split("_")[1],
                )
                for f in input.counts
            ]
        ).to_csv(output.csv, index=False)


rule ref_and_founder_nts:
    """Get nucleotide at each coding site for reference and clade founders."""
    input:
        coding_sites=rules.ref_coding_sites.output.csv,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
        fastas=lambda wc: [
            f"results/clade_founders_no_indels/{clade}.fa"
            for clade in clades_w_adequate_counts(wc)
        ],
    output:
        csv="results/coding_nts/coding_nts.csv",
    run:
        coding_sites=pd.read_csv(input.coding_sites)["site"].tolist()
        records = []
        for f_name in [input.ref_fasta, *input.fastas]:
            seq = str(Bio.SeqIO.read(f_name, "fasta").seq)
            name = os.path.splitext(os.path.basename(f_name))[0]
            for site in coding_sites:
                records.append((name, site, seq[site - 1]))
        os.makedirs(os.path.dirname(output.csv), exist_ok=True)
        pd.DataFrame.from_records(records, columns=["clade", "site", "nt"]).to_csv(
            output.csv, index=False,
        )


rule synonymous_mut_rates:
    """Compute and analyze rates and spectra of synonymous mutations."""
    input:
        csv=rules.aggregate_mutation_counts.output.csv,
        nb="notebooks/synonymous_mut_rates.ipynb",
    output:
        nb="results/synonymous_mut_rates/synonymous_mut_rates.ipynb",
        nb_html="results/synonymous_mut_rates/synonymous_mut_rates.html",
    params:
        synonymous_spectra_min_counts=config["synonymous_spectra_min_counts"],
        subset_order="{subset_order: " + str(list(config['sample_subsets'])) + "}",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p synonymous_spectra_min_counts {params.synonymous_spectra_min_counts} \
            -y "{params.subset_order}" \
            -p input_csv {input.csv}
        jupyter nbconvert {output.nb} --to html
        """
