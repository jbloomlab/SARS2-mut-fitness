# SARS-CoV-2 mutational spectrum

## Overview
This repository analyzes the mutational spectrum of human SARS-CoV-2 to understand shifts in the relative rates of different mutations over time.

This analysis was originally based off starting to replicate analyses [Neher (2022)](https://www.biorxiv.org/content/10.1101/2022.08.22.504731v1.full), but then checking the assumption that underlying synonymous mutation spectra were conserved across clades.

Study by Jesse Bloom, Richard Neher, Kelley Harris.

## Structure of repository and running the analysis
The analysis is entirely reproducible.

First build the [conda](https://docs.conda.io/) environment in [environment.yml](environment.yml).
This requires you to install [conda](https://docs.conda.io/), and then run:

    conda env create -f environment.yml

That command will create a `conda` environment named `SARS2-mut-spectrum` which you can activate with:

    conda activate SARS2-mut-spectrum

Then run the [snakemake](https://snakemake.readthedocs.io/) pipeline in [Snakefile](Snakefile), which reads its configuration from [config.yaml](config.yaml) by running:

    snakemake -j <n_cpus>

where `<n_cpus>` is the number of CPUs to use.

Note that the pipeline uses Python scripts in [./scripts/](scripts) and Jupyter notebooks in [./notebooks/](notebooks).

The created files are placed in [./results/](results).
Only some of those results files are tracked in this repo (others are too large to track).

## Overview of analysis
The analysis uses the pre-built mutation-annotated tree that is made available for use with the [UShER](https://usher-wiki.readthedocs.io/) package.
This tree contains all public access SARS-CoV-2 sequences, with mutations annotated on branches.
For the specific version of the tree used here, see the [config.yaml](config.yaml) file.

We analyze mutations grouping sequences at the level of [Nextstrain clades](https://clades.nextstrain.org/), which are already annotated on the pre-built mutation-annotated tree.
For each Nextstrain clade, we use the clade founder genotype manually defined by [Neher (2022)](https://www.biorxiv.org/content/10.1101/2022.08.22.504731v1.full) and available at the URL indicated in [config.yaml](config.yaml).

We perform some crucial filtering to remove spurious mutations as can arise from bad sequencing, base calling to reference, etc:

 1. We ignore all mutations on any branches with higher numbers of total mutations or mutations to either the reference or the founder of the clade in question. The exact settings for this filtering are in [config.yaml](config.yaml).

 2. We also ignore any branches with multiple mutations at the same codon.

 3. If a switch is set in [config.yaml](config.yaml), we ignore any mutations that are reversions from the clade founder to the reference, and also the reverse complement of these mutations. This is designed to remove missing bases called to reference, and also complements of those mutations induced by spurious nodes with such miscalls on downstream branches in the tree.

 4. We ignore all mutations at error-prone or problematic sites as specified in [config.yaml](config.yaml).

 5. We ignore any clades with small numbers of sequence samples as indicated in [config.yaml](config.yaml).

 6. Note also that indels are ignored, as they are not captured in the mutation-annotated tree.

We then analyze the mutation spectra for each clade.
For this analysis, we only consider synonymous mutations.
We also perform analyses for subsets of sequences from different regions (as specified in [config.yaml](config.yaml)) as well as for the the genome partitioned into halves--these analyses are designed to check that results are not due to artifacts related to sequencing pipelines or hotspots in the genome.
For all of these analyses, we only include subsets/partitions with at least the minimum number of mutations indicated in [config.yaml](config.yaml).

The pipeline in [Snakefile](Snakefile) does the various processing and filtering steps, and all of the final analysis is done by a single Jupyter notebook.
You can look at the HTML rendering of that Jupyter notebook at [results/synonymous_mut_rates/synonymous_mut_rates.html](results/synonymous_mut_rates/synonymous_mut_rates.html): download and open that notebook to look at the interactive `altair` plots.
