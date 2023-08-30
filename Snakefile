"""Top-level ``snakemake`` file that runs pipeline."""


import os
import textwrap

import pandas as pd

import yaml


configfile: "config.yaml"


with open(config["docs_plot_annotations"]) as f:
    docs_plot_annotations = yaml.safe_load(f)

# paths with results subdirectory for key tracked output files
results_files = [
    "synonymous_mut_rates/rates_by_clade.csv",
    "mutation_counts/aggregated.csv",
    "expected_vs_actual_mut_counts/expected_vs_actual_mut_counts.csv",
    "clade_founder_nts/clade_founder_nts.csv",
    "clade_founder_nts/clade_founder_aas.csv",
    "aa_fitness/aamut_fitness_all.csv",
    "aa_fitness/aamut_fitness_by_clade.csv",
    "aa_fitness/aamut_fitness_by_subset.csv",
    "aa_fitness/aa_fitness.csv",
    "aa_fitness/aa_fitness.json",
    "aa_fitness/aa_fitness.json.gz",
    "nt_fitness/ntmut_fitness_all.csv",
    "nt_fitness/ntmut_fitness_by_clade.csv",
    "nt_fitness/ntmut_fitness_by_subset.csv",
    "nt_fitness/nt_fitness.csv",
    "nt_fitness/synonymous_constraint_figure.pdf",
    "comparator_studies/comparator_corr.html",
    "dms-viz/mut_fitness.json",
    *[f"dms/{dms_dataset}/processed.csv" for dms_dataset in config["dms_datasets"]],
]


rule all:
    """Target rule with desired output files."""
    input:
        expand(
            [os.path.join("results_{mat}", f) for f in results_files],
            mat=config["mat_trees"],
        ),
        expand(
            [os.path.join("results", f) for f in results_files],
            mat=config["mat_trees"],
        ),
        expand(
            "docs/{mat}/{plot}.html",
            plot=list(docs_plot_annotations["plots"]) + ["index"],
            mat=config["mat_trees"],
        ),
        expand(
            "docs/{plot}.html",
            plot=list(docs_plot_annotations["plots"]) + ["index"],
            mat=config["mat_trees"],
        ),


rule get_mat_tree:
    """Get the pre-built mutation-annotated tree."""
    params:
        url=lambda wc: config["mat_trees"][wc.mat],
    output:
        mat="results_{mat}/mat/mat_tree.pb.gz",
    shell:
        "curl {params.url} > {output.mat}"


rule get_ref_fasta:
    """Get the reference FASTA."""
    params:
        url=config["ref_fasta"],
    output:
        ref_fasta="results_{mat}/ref/ref.fa",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_fasta}"


rule get_ref_gtf:
    """Get the reference GTF."""
    params:
        url=config["ref_gtf"],
    output:
        ref_gtf="results_{mat}/ref/original_ref.gtf",
    shell:
        "wget -O - {params.url} | gunzip -c > {output.ref_gtf}"


rule edit_ref_gtf:
    """Edit the reference GTF with manual additions."""
    input:
        gtf=rules.get_ref_gtf.output.ref_gtf,
    output:
        gtf="results_{mat}/ref/edited_ref.gtf",
    params:
        edits=config["add_to_ref_gtf"],
    notebook:
        "notebooks/edit_ref_gtf.py.ipynb"


rule ref_coding_sites:
    """Get all sites in reference that are part of a coding sequence."""
    input:
        gtf=rules.edit_ref_gtf.output.gtf,
        fasta=rules.get_ref_fasta.output.ref_fasta,
    output:
        csv="results_{mat}/ref/coding_sites.csv",
    script:
        "scripts/ref_coding_sites.py"


checkpoint mat_samples:
    """Get all samples in mutation-annotated tree with their dates and clades."""
    input:
        mat=rules.get_mat_tree.output.mat,
    output:
        csv="results_{mat}/mat/samples.csv",
        clade_counts="results_{mat}/mat/sample_clade_counts.csv",
    params:
        min_clade_samples=config["min_clade_samples"],
    script:
        "scripts/mat_samples.py"


def clades_w_adequate_counts(wc):
    """Return list of all clades with adequate sample counts."""
    return [
        clade
        for clade in (
            pd.read_csv(checkpoints.mat_samples.get(**wc).output.clade_counts)
            .query("adequate_sample_counts")["nextstrain_clade"]
            .tolist()
        )
        if clade not in config["clades_to_exclude"]
    ]


rule samples_by_clade_subset:
    """Get samples in mutation-annotated tree by nextstrain clade and subset."""
    input:
        csv=rules.mat_samples.output.csv,
    output:
        txt="results_{mat}/mat_by_clade_subset/{clade}_{subset}.txt",
    params:
        match_regex=lambda wc: config["sample_subsets"][wc.subset],
    run:
        (
            pd.read_csv(input.csv)
            .query("nextstrain_clade == @wildcards.clade")
            .query(f"sample.str.match('{params.match_regex}')")["sample"]
            .to_csv(output.txt, index=False, header=False)
        )


rule mat_clade_subset:
    """Get mutation-annotated tree for just a clade and subset."""
    input:
        mat=rules.get_mat_tree.output.mat,
        samples=rules.samples_by_clade_subset.output.txt,
    output:
        mat="results_{mat}/mat_by_clade_subset/{clade}_{subset}.pb",
    shell:
        """
        if [ -s {input.samples} ]; then
            echo "Extracting samples from {input.samples}"
            matUtils extract -i {input.mat} -s {input.samples} -O -o {output.mat}
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
        ref_gtf=rules.edit_ref_gtf.output.gtf,
    output:
        tsv="results_{mat}/mat_by_clade_subset/{clade}_{subset}_mutations.tsv",
    shell:
        """
        matUtils summary \
            -i {input.mat} \
            -g {input.ref_gtf} \
            -f {input.ref_fasta} \
            -t {output.tsv}
        """


rule mat_sample_path:
    """Get sample paths on MAT for clade. Use for nucleotide mutations."""
    # for explanation of why we need this rule in addition to `translate_mat`:
    # https://github.com/yatisht/usher/issues/336#issuecomment-1490764515
    input:
        mat=rules.mat_clade_subset.output.mat,
    output:
        tsv=temp("results_{mat}/mat_by_clade_subset/{clade}_{subset}_sample_paths.tsv"),
    shell:
        "matUtils extract -i {input.mat} --sample-paths {output.tsv}"


rule sample_path_to_nt_mutations:
    """Get all nucleotide mutations on tree from sample paths."""
    # for explanation of why we need this rule in addition to `translate_mat`:
    # https://github.com/yatisht/usher/issues/336#issuecomment-1490764515
    input:
        tsv=rules.mat_sample_path.output.tsv,
    output:
        csv="results_{mat}/mat_by_clade_subset/{clade}_{subset}_nt_mutations.csv",
    script:
        "scripts/sample_path_to_nt_mutations.py"


rule clade_founder_jsons:
    """Get JSONs with nexstrain clade founders (indels not included)."""
    params:
        **config["clade_founders"],
    output:
        neher_json="results_{mat}/clade_founders_no_indels/clade_founders_neher.json",
        roemer_json="results_{mat}/clade_founders_no_indels/clade_founders_roemer.json",
    shell:
        """
        curl {params.neher_json} > {output.neher_json}
        curl {params.roemer_json} > {output.roemer_json}
        """


rule clade_founder_fasta_and_muts:
    """Get FASTA and mutations for nextstrain clade founder (indels not included)."""
    input:
        neher_json=rules.clade_founder_jsons.output.neher_json,
        roemer_json=rules.clade_founder_jsons.output.roemer_json,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
    params:
        roemer_nextstrain_to_pango=config["roemer_nextstrain_to_pango"],
    output:
        fasta="results_{mat}/clade_founders_no_indels/{clade}.fa",
        muts="results_{mat}/clade_founders_no_indels/{clade}_ref_to_founder_muts.csv",
    script:
        "scripts/clade_founder_fasta.py"


rule site_mask_vcf:
    """Get the site mask VCF."""
    output:
        vcf="results_{mat}/site_mask/site_mask.vcf",
    params:
        url=config["site_mask_vcf"],
    shell:
        "curl {params.url} > {output.vcf}"


rule site_mask:
    """Convert site mask VCF to CSV."""
    input:
        vcf=rules.site_mask_vcf.output.vcf,
    output:
        csv="results_{mat}/site_mask/site_mask.csv",
    script:
        "scripts/site_mask.py"


rule count_mutations:
    """Count mutations, excluding branches with too many mutations or reversions."""
    input:
        tsv=rules.translate_mat.output.tsv,
        nt_mut_csv=rules.sample_path_to_nt_mutations.output.csv,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
        clade_founder_fasta=rules.clade_founder_fasta_and_muts.output.fasta,
        ref_to_founder_muts=rules.clade_founder_fasta_and_muts.output.muts,
        usher_masked_sites=config["usher_masked_sites"],
        site_mask=rules.site_mask.output.csv,
    output:
        csv="results_{mat}/mutation_counts/{clade}_{subset}.csv",
    params:
        max_nt_mutations=config["max_nt_mutations"],
        max_reversions_to_ref=config["max_reversions_to_ref"],
        max_reversions_to_clade_founder=config["max_reversions_to_clade_founder"],
        exclude_ref_to_founder_muts=config["exclude_ref_to_founder_muts"],
        sites_to_exclude=config["sites_to_exclude"],
        site_include_range=config["site_include_range"],
    log:
        notebook="results_{mat}/mutation_counts/{clade}_{subset}_count_mutations.ipynb",
    notebook:
        "notebooks/count_mutations.py.ipynb"


rule clade_founder_nts:
    """Get nucleotide at each coding site for clade founders."""
    input:
        coding_sites=rules.ref_coding_sites.output.csv,
        fastas=lambda wc: [
            f"results_{wc.mat}/clade_founders_no_indels/{clade}.fa"
            for clade in clades_w_adequate_counts(wc)
        ],
    output:
        csv="results_{mat}/clade_founder_nts/clade_founder_nts.csv",
    script:
        "scripts/clade_founder_nts.py"


rule aggregate_mutation_counts:
    """Aggregate the mutation counts for all clades and subsets."""
    input:
        clade_founder_nts=rules.clade_founder_nts.output.csv,
        counts=lambda wc: [
            f"results_{wc.mat}/mutation_counts/{clade}_{subset}.csv"
            for clade in clades_w_adequate_counts(wc)
            for subset in config["sample_subsets"]
        ],
    output:
        csv="results_{mat}/mutation_counts/aggregated.csv",
    script:
        "scripts/aggregate_mutation_counts.py"


rule synonymous_mut_rates:
    """Compute and analyze rates and spectra of synonymous mutations."""
    input:
        mutation_counts_csv=rules.aggregate_mutation_counts.output.csv,
        clade_founder_nts_csv=rules.clade_founder_nts.output.csv,
    output:
        rates_by_clade="results_{mat}/synonymous_mut_rates/rates_by_clade.csv",
        rates_plot="results_{mat}/synonymous_mut_rates/mut_rates.html",
    params:
        synonymous_spectra_min_counts=config["synonymous_spectra_min_counts"],
        sample_subsets=config["sample_subsets"],
        clade_synonyms=config["clade_synonyms"],
    log:
        notebook="results_{mat}/synonymous_mut_rates/synonymous_mut_rates.ipynb",
    notebook:
        "notebooks/synonymous_mut_rates.ipynb"


rule expected_mut_counts:
    """Compute expected mutation counts from synonymous mutation rates and counts."""
    input:
        rates_by_clade=rules.synonymous_mut_rates.output.rates_by_clade,
        clade_founder_nts_csv=rules.clade_founder_nts.output.csv,
    output:
        expected_counts="results_{mat}/expected_mut_counts/expected_mut_counts.csv",
    log:
        notebook="results_{mat}/expected_mut_counts/expected_mut_counts.ipynb",
    notebook:
        "notebooks/expected_mut_counts.ipynb"


rule aggregate_mutations_to_exclude:
    """Aggregate the set of all mutations to exclude for each clade."""
    input:
        muts_to_exclude=lambda wc: [
            f"results_{wc.mat}/clade_founders_no_indels/{clade}_ref_to_founder_muts.csv"
            for clade in clades_w_adequate_counts(wc)
        ],
        usher_masked_sites=config["usher_masked_sites"],
        site_mask=rules.site_mask.output.csv,
        ref_fasta=rules.get_ref_fasta.output.ref_fasta,
    output:
        csv="results_{mat}/expected_vs_actual_mut_counts/mutations_to_exclude.csv",
    params:
        clades=lambda wc: clades_w_adequate_counts(wc),
        sites_to_exclude=config["sites_to_exclude"],
        site_include_range=config["site_include_range"],
        exclude_ref_to_founder_muts=config["exclude_ref_to_founder_muts"],
    script:
        "scripts/aggregate_mutations_to_exclude.py"


rule merge_expected_and_actual_counts:
    """Merge expected and actual counts."""
    input:
        expected=rules.expected_mut_counts.output.expected_counts,
        actual=rules.aggregate_mutation_counts.output.csv,
        muts_to_exclude=rules.aggregate_mutations_to_exclude.output.csv,
    output:
        csv="results_{mat}/expected_vs_actual_mut_counts/expected_vs_actual_mut_counts.csv",
    log:
        notebook="results_{mat}/expected_vs_actual_mut_counts/merge_expected_and_actual_counts.ipynb",
    notebook:
        "notebooks/merge_expected_and_actual_counts.py.ipynb"


rule summarize_expected_vs_actual:
    """Summarize expected vs actual across mutations."""
    input:
        csv=rules.merge_expected_and_actual_counts.output.csv,
    output:
        chart="results_{mat}/expected_vs_actual_mut_counts/avg_counts.html",
    log:
        notebook="results_{mat}/expected_vs_actual_mut_counts/summarize_expected_vs_actual.ipynb",
    notebook:
        "notebooks/summarize_expected_vs_actual.py.ipynb"


rule aamut_fitness:
    """Fitness effects from expected vs actual counts for amino-acid mutations."""
    input:
        csv=rules.merge_expected_and_actual_counts.output.csv,
    output:
        aamut_all="results_{mat}/aa_fitness/aamut_fitness_all.csv",
        aamut_by_clade="results_{mat}/aa_fitness/aamut_fitness_by_clade.csv",
        aamut_by_subset="results_{mat}/aa_fitness/aamut_fitness_by_subset.csv",
    params:
        orf1ab_to_nsps=config["orf1ab_to_nsps"],
        fitness_pseudocount=config["fitness_pseudocount"],
        gene_overlaps=config["gene_overlaps"],
    log:
        notebook="results_{mat}/aa_fitness/aamut_fitness.ipynb",
    notebook:
        "notebooks/aamut_fitness.py.ipynb"


rule ntmut_fitness:
    """Fitness effects from expected vs actual counts for nucleotide mutations."""
    input:
        csv=rules.merge_expected_and_actual_counts.output.csv,
    output:
        ntmut_all="results_{mat}/nt_fitness/ntmut_fitness_all.csv",
        ntmut_by_clade="results_{mat}/nt_fitness/ntmut_fitness_by_clade.csv",
        ntmut_by_subset="results_{mat}/nt_fitness/ntmut_fitness_by_subset.csv",
    params:
        fitness_pseudocount=config["fitness_pseudocount"],
    log:
        notebook="results_{mat}/nt_fitness/ntmut_fitness.ipynb",
    notebook:
        "notebooks/ntmut_fitness.py.ipynb"


rule aa_fitness:
    """Fitnesses of different amino acids across clades."""
    input:
        aamut_fitness=rules.aamut_fitness.output.aamut_all,
    output:
        aa_fitness="results_{mat}/aa_fitness/aa_fitness.csv",
    log:
        notebook="results_{mat}/aa_fitness/aa_fitness.ipynb",
    notebook:
        "notebooks/aa_fitness.py.ipynb"


rule nt_fitness:
    """Fitnesses of different nucleotides across clades."""
    input:
        ntmut_fitness=rules.ntmut_fitness.output.ntmut_all,
    output:
        nt_fitness="results_{mat}/nt_fitness/nt_fitness.csv",
    log:
        notebook="results_{mat}/nt_fitness/nt_fitness.ipynb",
    notebook:
        "notebooks/nt_fitness.py.ipynb"


rule clade_founder_aas:
    """Get clade-founder amino acids."""
    input:
        clade_founder_nts=rules.clade_founder_nts.output.csv,
    params:
        orf1ab_to_nsps=config["orf1ab_to_nsps"],
        clade_synonyms=config["clade_synonyms"],
    output:
        clade_founder_aas="results_{mat}/clade_founder_nts/clade_founder_aas.csv",
    notebook:
        "notebooks/clade_founder_aas.py.ipynb"


rule analyze_aa_fitness:
    """Analyze and plot amino-acid mutation fitnesses."""
    input:
        aamut_all=rules.aamut_fitness.output.aamut_all,
        aamut_by_subset=rules.aamut_fitness.output.aamut_by_subset,
        aamut_by_clade=rules.aamut_fitness.output.aamut_by_clade,
        aafitness=rules.aa_fitness.output.aa_fitness,
        clade_founder_nts=rules.clade_founder_nts.output.csv,
        ref_coding_sites=rules.ref_coding_sites.output.csv,
    params:
        min_expected_count=config["min_expected_count"],
        clade_corr_min_count=config["clade_corr_min_count"],
        init_ref_clade=config["aa_fitness_init_ref_clade"],
        clade_synonyms=config["clade_synonyms"],
        heatmap_minimal_domain=config["aa_fitness_heatmap_minimal_domain"],
        orf1ab_to_nsps=config["orf1ab_to_nsps"],
    output:
        outdir=directory("results_{mat}/aa_fitness/plots"),
    log:
        notebook="results_{mat}/aa_fitness/analyze_aa_fitness.ipynb",
    notebook:
        "notebooks/analyze_aa_fitness.py.ipynb"


rule get_dms_dataset:
    """Get a deep mutational scanning dataset."""
    params:
        url=lambda wc: config["dms_datasets"][wc.dms_dataset]["url"],
    output:
        raw_data="results_{mat}/dms/{dms_dataset}/raw.csv",
    shell:
        "curl {params.url} > {output.raw_data}"


rule process_dms_dataset:
    """Process a deep mutational scanning dataset to fitness estimates."""
    input:
        unpack(
            lambda wc: (
                {"wt_seq": config["dms_datasets"][wc.dms_dataset]["wt_seq"]}
                if "wt_seq" in config["dms_datasets"][wc.dms_dataset]
                else {}
            )
        ),
        raw_data=rules.get_dms_dataset.output.raw_data,
    output:
        processed="results_{mat}/dms/{dms_dataset}/processed.csv",
    log:
        notebook="results_{mat}/dms/{dms_dataset}/process_{dms_dataset}.ipynb",
    notebook:
        "notebooks/process_{wildcards.dms_dataset}.ipynb"


rule fitness_dms_corr:
    """Correlate the fitness estimates with those from deep mutational scanning."""
    input:
        **{
            dms_dataset: os.path.join(
                "results_{mat}", "dms", dms_dataset, "processed.csv"
            )
            for dms_dataset in config["dms_datasets"]
        },
        aafitness=rules.aa_fitness.output.aa_fitness,
        neher_fitness=config["neher_fitness"],
    output:
        plotsdir=directory("results_{mat}/fitness_dms_corr/plots"),
    params:
        min_expected_count=config["min_expected_count"],
        dms_datasets=config["dms_datasets"],
    log:
        notebook="results_{mat}/fitness_dms_corr/fitness_dms_corr.ipynb",
    notebook:
        "notebooks/fitness_dms_corr.py.ipynb"


rule clade_fixed_muts:
    """Analyze mutations fixed in each clade."""
    input:
        aafitness=rules.aa_fitness.output.aa_fitness,
        aamut_by_clade=rules.aamut_fitness.output.aamut_by_clade,
        clade_founder_nts_csv=rules.clade_founder_nts.output.csv,
    output:
        fixed_muts_chart="results_{mat}/clade_fixed_muts/clade_fixed_muts.html",
        fixed_muts_hist="results_{mat}/clade_fixed_muts/clade_fixed_muts_hist.html",
    params:
        min_expected_count=config["min_expected_count"],
        ref=config["clade_fixed_muts_ref"],
        orf1ab_to_nsps=config["orf1ab_to_nsps"],
    log:
        notebook="results_{mat}/clade_fixed_muts/clade_fixed_muts.ipynb",
    notebook:
        "notebooks/clade_fixed_muts.py.ipynb"


rule fitness_vs_terminal:
    """Analyze fitness effects of mutations vs terminal / non-terminal node counts."""
    input:
        aamut_all=rules.aamut_fitness.output.aamut_all,
    output:
        chart="results_{mat}/fitness_vs_terminal/fitness_vs_terminal.html",
    params:
        min_expected_count=config["min_expected_count"],
        min_actual_count=config["terminal_min_actual_count"],
        pseudocount=config["terminal_pseudocount"],
    log:
        notebook="results_{mat}/fitness_vs_terminal/fitness_vs_terminal.ipynb",
    notebook:
        "notebooks/fitness_vs_terminal.py.ipynb"


rule synonymous_figures:
    """Plot of synonymous selection."""
    input:
        fitness=rules.ntmut_fitness.output.ntmut_all,
    output:
        synonymous_figure="results_{mat}/nt_fitness/synonymous_constraint_figure.pdf",
    params:
        min_expected_count=config["min_expected_count"],
    shell:
        """
        python scripts/noncoding_constraints.py \
            --fitness {input.fitness} \
            --output {output.synonymous_figure} \
            --min_expected_count {params.min_expected_count}
        """


rule correlate_mats:
    """Correlate mutation effects for different MATs (sequence sets)."""
    input:
        aa_fitnesses=expand(rules.aa_fitness.output.aa_fitness, mat=config["mat_trees"]),
    output:
        fitness_corrs_chart="results_{mat}/mat_corrs/mat_aa_fitness_correlations.html",
    params:
        mats=list(config["mat_trees"]),
        min_expected_count=config["min_expected_count"],
    log:
        notebook="results_{mat}/mat_corrs/correlate_mats.ipynb",
    notebook:
        "notebooks/correlate_mats.py.ipynb"


rule get_dnds_data:
    """Get data for dN/dS analysis."""
    params:
        url=config["dnds"],
    output:
        csv="results_{mat}/dnds/dnds_data.csv",
    shell:
        "curl {params.url} > {output.csv}"


rule analyze_dnds:
    """Analyze dN/dS versus amino-acid fitness and DMS data."""
    input:
        **{
            dms_dataset: os.path.join(
                "results_{mat}", "dms", dms_dataset, "processed.csv"
            )
            for dms_dataset in config["dms_datasets"]
        },
        dnds=rules.get_dnds_data.output.csv,
        aa_fitness=rules.aa_fitness.output.aa_fitness,
    output:
        corr_html="results_{mat}/dnds/dnds_corr.html",
    params:
        min_expected_count=config["min_expected_count"],
        dms_datasets=config["dms_datasets"],
    log:
        notebook="results_{mat}/dnds/analyze_dnds.ipynb",
    notebook:
        "notebooks/analyze_dnds.py.ipynb"


rule analyze_comparator_studies:
    """Analyze comparator studies versus fitness estimates and DMS data."""
    input:
        **{
            study: f"data/comparator_studies/{study}.csv"
            for study in config["comparator_studies"]
        },
        **{
            dms_dataset: os.path.join(
                "results_{mat}", "dms", dms_dataset, "processed.csv"
            )
            for dms_dataset in config["dms_datasets"]
        },
        aa_fitness=rules.aa_fitness.output.aa_fitness,
    output:
        corr_html="results_{mat}/comparator_studies/comparator_corr.html",
    params:
        min_expected_count=config["min_expected_count"],
        dms_datasets=config["dms_datasets"],
        comparator_studies=config["comparator_studies"],
    log:
        notebook="results_{mat}/comparator_studies/analyze_comparator_studies.ipynb",
    notebook:
        "notebooks/analyze_comparator_studies.py.ipynb"


rule aggregate_plots_for_docs:
    """Aggregate plots to include in GitHub pages docs."""
    input:
        aa_fitness_plots_dir=rules.analyze_aa_fitness.output.outdir,
        dms_corr_plotsdir=rules.fitness_dms_corr.output.plotsdir,
        rates_plot=rules.synonymous_mut_rates.output.rates_plot,
        clade_fixed_muts=rules.clade_fixed_muts.output.fixed_muts_chart,
        clade_fixed_hist=rules.clade_fixed_muts.output.fixed_muts_hist,
        fitness_vs_terminal=rules.fitness_vs_terminal.output.chart,
        avg_counts=rules.summarize_expected_vs_actual.output.chart,
        mat_corrs=rules.correlate_mats.output.fitness_corrs_chart,
        dnds_corr=rules.analyze_dnds.output.corr_html,
        comparator_corr=rules.analyze_comparator_studies.output.corr_html,
    output:
        expand(
            "results_{{mat}}/plots_for_docs/{plot}.html",
            plot=docs_plot_annotations["plots"],
        ),
    params:
        plotsdir="results_{mat}/plots_for_docs",
    shell:
        """
        mkdir -p {params.plotsdir}
        rm -f {params.plotsdir}/*
        cp {input.aa_fitness_plots_dir}/*.html {params.plotsdir}
        cp {input.dms_corr_plotsdir}/*.html {params.plotsdir}
        cp {input.rates_plot} {params.plotsdir}
        cp {input.clade_fixed_muts} {params.plotsdir}
        cp {input.clade_fixed_hist} {params.plotsdir}
        cp {input.fitness_vs_terminal} {params.plotsdir}
        cp {input.avg_counts} {params.plotsdir}
        cp {input.mat_corrs} {params.plotsdir}
        cp {input.dnds_corr} {params.plotsdir}
        cp {input.comparator_corr} {params.plotsdir}
        """


rule format_plot_for_docs:
    """Format a specific plot for the GitHub pages docs."""
    input:
        plot=os.path.join(rules.aggregate_plots_for_docs.params.plotsdir, "{plot}.html"),
        script="scripts/format_altair_html.py",
    output:
        plot="docs/{mat}/{plot}.html",
        markdown=temp("results_{mat}/plots_for_docs/{plot}.md"),
    params:
        annotations=lambda wc: docs_plot_annotations["plots"][wc.plot],
        url=config["docs_url"],
        legend_suffix=docs_plot_annotations["legend_suffix"],
    shell:
        """
        echo "## {params.annotations[title]}\n" > {output.markdown}
        echo "{params.annotations[legend]}\n\n" >> {output.markdown}
        echo "{params.legend_suffix}" >> {output.markdown}
        echo "\nThis plot is for the {wildcards.mat} dataset." >> {output.markdown}
        echo "[Here](index.html) are all plots for that dataset." >> {output.markdown}
        python {input.script} \
            --chart {input.plot} \
            --markdown {output.markdown} \
            --site {params.url} \
            --title "{params.annotations[title]}" \
            --description "{params.annotations[title]}" \
            --output {output.plot}
        """


rule docs_index:
    """Write index for GitHub Pages docs for each MAT."""
    output:
        html="docs/{mat}/index.html",
    params:
        plot_annotations=docs_plot_annotations,
        mat_trees=list(config["mat_trees"]),
    script:
        "scripts/docs_index.py"


rule current_docs_index:
    """Write index for GitHub Pages docs for current MAT."""
    output:
        html="docs/index.html",
    params:
        plot_annotations=docs_plot_annotations,
        mat_trees=list(config["mat_trees"]),
        current_mat=config["current_mat"],
    script:
        "scripts/docs_index.py"


rule cp_current_mat_results:
    """Copy the current MAT results to `./results/`."""
    input:
        [os.path.join(f"results_{config['current_mat']}", f) for f in results_files],
    output:
        [os.path.join("results", f) for f in results_files],
        subdir=directory("results"),
    params:
        input_subdir=f"results_{config['current_mat']}",
    shell:
        """
        rm -rf {output.subdir}
        cp -r {params.input_subdir} {output.subdir}
        """


rule cp_current_mat_docs:
    """Copy the current MAT docs to `./docs/` for default GitHub pages display."""
    input:
        expand(
            os.path.join("docs", config["current_mat"], "{plot}.html"),
            plot=list(docs_plot_annotations["plots"]),
        ),
    output:
        expand(
            os.path.join("docs", "{plot}.html"),
            plot=list(docs_plot_annotations["plots"]),
        ),
    shell:
        "cp {input} docs"


rule export_fitness_to_json:
    """Export amino-acid fitness to JSON with metadata"""
    input:
        aa_fitness=rules.aa_fitness.output.aa_fitness,
        clade_founder_aas=rules.clade_founder_aas.output.clade_founder_aas,
    output:
        aa_fitness_json="results_{mat}/aa_fitness/aa_fitness.json",
        aa_fitness_json_gz="results_{mat}/aa_fitness/aa_fitness.json.gz",
    params:
        min_expected_count=config["min_expected_count"],
        citation=config["citation"],
        authors=config["authors"],
        source=config["source"],
        description=config["description"],
    script:
        "scripts/export_fitness_to_json.py"


# Format the data for each protein
rule format_fitness_for_dms_viz:
    input:
        aa_fitness=rules.aa_fitness.output.aa_fitness,
        clade_founder_aas=rules.clade_founder_aas.output.clade_founder_aas,
        structure_info="data/proteins.csv",
    output:
        fitness_df="results_{mat}/dms-viz/{protein}/{protein}_fitness.csv",
        sitemap_df="results_{mat}/dms-viz/{protein}/{protein}_sitemap.csv",
    script:
        "scripts/format-data-for-dms-viz.py"


# Create JSON files for each protein
rule create_dms_viz_json:
    input:
        fitness_df=("results_{mat}/dms-viz/{protein}/{protein}_fitness.csv"),
        sitemap_df=("results_{mat}/dms-viz/{protein}/{protein}_sitemap.csv"),
    output:
        os.path.join("results_{mat}/dms-viz/", "{protein}", "{protein}.json"),
    params:
        name=lambda wildcards: wildcards.protein,
        structure=lambda wildcards: proteins.loc[
            proteins["selection"] == wildcards.protein, "pdb"
        ].item(),
        include_chains=lambda wildcards: proteins.loc[
            proteins["selection"] == wildcards.protein, "dataChains"
        ].item(),
        exclude_chains=lambda wildcards: proteins.loc[
            proteins["selection"] == wildcards.protein, "excludedChains"
        ].item(),
        description=lambda wildcards: proteins.loc[
            proteins["selection"] == wildcards.protein, "description"
        ].item(),
        title=lambda wildcards: proteins.loc[
            proteins["selection"] == wildcards.protein, "title"
        ].item(),
        filter_cols={"expected_count": "Expected Count"},
        filter_limits={"expected_count": [0, 100]},
        tooltip_cols={"expected_count": "Expected Count"},
        metric="fitness",
        metric_name="Fitness",
    conda: "envs/configure-dms-viz.yml"
    shell:
        """
        configure-dms-viz \
            --input {input.fitness_df} \
            --name "{params.name}" \
            --sitemap {input.sitemap_df} \
            --metric {params.metric} \
            --structure {params.structure} \
            --metric-name {params.metric_name} \
            --output {output} \
            --alphabet "RKHDEQNSTYWFAILMVGPC*" \
            --included-chains "{params.include_chains}" \
            --excluded-chains "{params.exclude_chains}" \
            --filter-cols "{params.filter_cols}" \
            --filter-limits "{params.filter_limits}" \
            --tooltip-cols "{params.tooltip_cols}" \
            --exclude-amino-acids "*" \
            --description "{params.description}" \
            --title "{params.title}"
        """


# Combine JSON files into one
proteins = pd.read_csv("data/proteins.csv")
rule combine_dms_viz_jsons:
    input:
        input_files=lambda wildcards: [
            os.path.join(
                f"results_{wildcards.mat}/dms-viz/", protein, f"{protein}.json"
            )
            for protein in proteins.selection.unique()
        ],
    output:
        output_file="results_{mat}/dms-viz/mut_fitness.json",
    run:
        combined_data = {}
        for input_file in input.input_files:
            with open(input_file) as f:
                data = json.load(f)
                combined_data.update(data)
        with open(output.output_file, "w") as f:
            json.dump(combined_data, f)
