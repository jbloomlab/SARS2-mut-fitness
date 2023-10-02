# Fitness effects of SARS-CoV-2 amino-acid mutations estimated from observed versus expected mutation counts

## Overview
This repository estimates the fitness effects of mutations to all SARS-CoV-2 proteins as described in [this paper](https://academic.oup.com/ve/article/9/2/vead055/7265011).
It does this by analyzing mutations in millions human SARS-CoV-2 sequences, and quantifying how the observed counts of each mutation compares to the expected counts from the underlying neutral mutation rate.

The analysis is by [Jesse Bloom](https://scholar.google.com/citations?user=S12x_eQAAAAJ&hl=en) and [Richard Neher](https://neherlab.org/), and makes use of the SARS-CoV-2 mutation-annotated tree provided by the [developers of UShER](https://usher-wiki.readthedocs.io/).

## References
For details, see the following references:
 - The main paper describing this work is [Bloom and Neher (2023)](https://academic.oup.com/ve/article/9/2/vead055/7265011).
 - Secondary references:
    - The approach builds on ideas initially described in [Neher (2022)](https://doi.org/10.1093/ve/veac113).
    - The estimation of the neutral mutation rate is done as described in [Bloom, Beichman, Neher, & Harris (2023)](https://academic.oup.com/mbe/article/40/4/msad085/7113660).
    - The mutational data are extracted from publicly available SARS-CoV-2 sequences using the UShER package described by [Turakhia et al (2021)](https://www.nature.com/articles/s41588-021-00862-7).

## Interactive plots of results
The easiest way to access the results is through a set of interactive plots available at [https://jbloomlab.github.io/SARS2-mut-fitness](https://jbloomlab.github.io/SARS2-mut-fitness). You can also visualize the fitness effects of mutations in the context of a 3D protein structure [here]() using the web-tool `dms-viz`.

## CSV files with numerical results
Here are links to files with major numerical results:
 - Final estimates of relative [fitnesses of each amino acid](results/aa_fitness/aa_fitness.csv). For most purposes, you probably want to use this file for the final "best" estimate of the amino-acid fitnesses. Note that estimates are more accurate for amino acids with larger values of *expected_count*. Sites in ORF1ab / nsp proteins are listed both with the ORF1ab and nsp numbering.
 - Estimates of the effects of amino-acid mutations [aggregated across all clades](results/aa_fitness/aamut_fitness_all.csv), for each [individual clade](results/aa_fitness/aamut_fitness_by_clade.csv), and for just [subsets of sequences from specific countries](results/aa_fitness/aamut_fitness_by_subset.csv).
 - Estimates of the relative neutral mutation rates for different types of mutations for each clade, as estimated at four-fold degenerate synonymous sites, are [here](results/synonymous_mut_rates/rates_by_clade.csv).
 - The counts of each mutation at each site in each clade are [here](results/mutation_counts/aggregated.csv); note this file also contains excluded sites.
 - The observed counts of each mutation versus the counts expected from the underlying neutral mutation rate are [here](results/expected_vs_actual_mut_counts/expected_vs_actual_mut_counts.csv); note that this file include excluded sites.
 - The nucleotide identities at each position in the "founder" sequence for each viral clade are [here](results/clade_founder_nts/clade_founder_nts.csv). This file also indicates which sites are four-fold degenerate.
 - The deep mutational scanning results processed and aggregated from published experimental studies are in [this subdirectory](results/dms/).

## Different dataset versions
You can run the pipeline for multiple mutation-annotated tree datasets.
Specify the dataset to show by default in the interactive plots and to store results in [./results/](results) for as the `current_mat` in [config.yaml](config.yaml).
Other datasets specifies under the `mat_trees` key in [config.yaml](config.yaml) will also be analyzed, and their results go in subdirectories called `./results_{mat}/` and their interactive plots are available via the GitHub pages interactive rendering at [https://jbloomlab.github.io/SARS2-mut-fitness](https://jbloomlab.github.io/SARS2-mut-fitness) in a list at the bottom of the page.

## Structure of repository and running the analysis
The analysis is entirely reproducible from the code provided in this GitHub repository.

First build the [conda](https://docs.conda.io/) environment in [environment.yml](environment.yml).
This requires you to install [conda](https://docs.conda.io/), and then run:

    conda env create -f environment.yml

That command will create a `conda` environment named `SARS2-mut-fitness` which you can activate with:

    conda activate SARS2-mut-fitness

Then run the [snakemake](https://snakemake.readthedocs.io/) pipeline in [Snakefile](Snakefile), which reads its configuration from [config.yaml](config.yaml) by running:

    snakemake -j <n_cpus> --use-conda

where `<n_cpus>` is the number of CPUs to use.

Note that the pipeline uses Python scripts in [./scripts/](scripts) and Jupyter notebooks in [./notebooks/](notebooks).

The created files are placed in [./results/](results).
Only some of those results files are tracked in this repo (others are too large or numerous to track).
The output interactive HTML [altair](https://github.com/altair-viz/altair) plots are placed in [./docs/](docs) where they are displayed via GitHub pages at [https://jbloomlab.github.io/SARS2-mut-fitness](https://jbloomlab.github.io/SARS2-mut-fitness)

In addition to the configuration in [config.yaml](config.yaml), there is also some input data / specifications in [./data/](data).
This includes the file[data/docs_plot_annotations.yaml usher_masked_sites.yaml](data/docs_plot_annotations.yaml usher_masked_sites.yaml) that specifies how to layout and label the plots rendered on GitHub pages.

## Analysis
Basic steps, performed for each Nextstrain clade:

 - Count occurrences of all mutations along the [UShER](https://usher-wiki.readthedocs.io/) tree with some filtering to remove spurious ones.

 - Estimate the underlying mutation rates using 4-fold degenerate synonymous sites.

 - Based on the mutation rates and total number of observed synonymous mutations at those sites, calculate **expected** number of occurrences of each nucleotide mutation.

 - Compare observed to expected mutation counts to estimate the fitness effects of amino-acid mutations.

 - Make overall estimates of amino-acid fitnesses by aggregating mutation-effect estimates across clades.

 - Some additional plotting and analyses, as well as comparison to experimental measurements from deep mutational scanning studies.

### Counting mutations
Mutation counts are extracted from the pre-built mutation-annotated tree that is made available for use with the [UShER](https://usher-wiki.readthedocs.io/) package.
This tree contains all public access SARS-CoV-2 sequences, with mutations annotated on branches.
For the specific version of the tree used here, see the [config.yaml](config.yaml) file.

We analyze mutations grouping sequences at the level of [Nextstrain clades](https://clades.nextstrain.org/), which are already annotated on the pre-built mutation-annotated tree.
For each Nextstrain clade, we use the clade founder genotype manually defined by [Neher (2022)](https://doi.org/10.1093/ve/veac113) and available at the URL indicated in [config.yaml](config.yaml), or for newer clades the sequences from Cornelius Roemer [here](https://github.com/corneliusroemer/pango-sequences/blob/main/data/pango-consensus-sequences_summary.json).

We perform some crucial filtering to remove spurious mutations as can arise from bad sequencing, base calling to reference, etc:

 1. We ignore all mutations on any branches with high numbers of total mutations, or mutations to either the reference or the founder sequence of the clade in question. The exact settings for this filtering are in [config.yaml](config.yaml).

 2. We ignore any branches with multiple mutations at the same codon.

 3. If a switch is set in [config.yaml](config.yaml) (it currently is), we specify to exclude any mutations that are reversions from the clade founder to the reference, and also the reverse complement of these mutations. This is designed to remove missing bases called to reference, and also complements of those mutations induced by spurious nodes with such miscalls on downstream branches in the tree.

 4. We specify to exclude all mutations at error-prone or problematic sites as manually specified in [config.yaml](config.yaml). We also exclude the masked sites in the `UShER` pre-built tree specified [here](https://github.com/W-L/ProblematicSites_SARS-CoV2/). We also ignore for specific clades any mutations at sites that are masked in the `UShER` pipeline for those clades as specified in [data/usher_masked_sites.yaml](data/usher_masked_sites.yaml).

 5. We ignore any clades with small numbers of sequence samples as indicated in [config.yaml](config.yaml) as these are expected to have too much noise.

 6. Note also that indels are ignored, as they are not captured in the mutation-annotated tree.

 7. Although the main analysis here uses the total counts of each mutation, we also keep track of how many of these counts are on non-terminal (interior) branches of the tree versus terminal (tip) branches, and also the mean log descendants defined as the log of the number of leaves sharing mutations from each mutated branch in the UShER tree (with log zero values set to one).

The above mutation counts both for all sequences for a clade, and for the sample subsets defined in [config.yaml](config.yaml) are stored in [results/mutation_counts/aggregated.csv](results/mutation_counts/aggregated.csv).
Note that mutations are annotated by the protein(s) they affect, if they are synonymous, if they are at 4-fold degenerate sites, if they are at an excluded site, etc.

### Analysis of 4-fold degenerate synonymous mutation spectrum / rates

We then analyze the mutation spectra and rate of mutations.
For this analysis, we only consider synonymous mutations at sites (third codon positions) that are four-fold degenerate in the founder sequence for each clade.
The file [results/clade_founder_nts/clade_founder_nts.csv](results/clade_founder_nts/clade_founder_nts.csv) specifies which these sites are for each clade.

Specifically, we determine the relative fraction of all 4-fold synonymous mutations that are each type of nucleotide change, and also the relative **rates** of the different types of mutations, which are just computed as the fraction of all mutations of that type normalized by the composition of the sequence at these 4-fold degenerate sites.
These rates for each clade with sufficient counts are written to the file [results/synonymous_mut_rates/rates_by_clade.csv](results/synonymous_mut_rates/rates_by_clade.csv).

We also perform analyses for subsets of sequences from different regions (as specified in [config.yaml](config.yaml)) as well as for the the genome partitioned into halves--these analyses are designed to check that results are not due to artifacts related to sequencing pipelines or hotspots in the genome.
For all of these analyses, we only include subsets/partitions with at least the minimum number of mutations indicated in [config.yaml](config.yaml).

Most of the analysis of the synonymous mutation spectrum is done by [notebooks/synonymous_mut_rates.ipynb](notebooks/synonymous_mut_rates.ipynb).

A plot of the neutral mutation rates is available at [https://jbloomlab.github.io/SARS2-mut-fitness](https://jbloomlab.github.io/SARS2-mut-fitness)

### Computation of "expected" number of occurrences for each mutation

We next compute the "expected" number of observations of each mutation in the absence of selection based on the underlying mutation rates.
Specifically, above we have computed the relative rates of each type of mutation at 4-fold degenerate sites (they are normalized to the frequency of the parent nucleotide).
To get the expected numbers, for each clade, we first compute $T$ that satisfies 
$$N_s = T \sum\limits_{nt_1} s_{nt_1} \sum\limits_{nt_2 \ne nt_1} r_{nt_1\rightarrow nt_2}$$
where $N_s$ is the total number of mutations at 4-fold degenerate synonymous sites observed for the clade, $s_{nt}$ is the number of 4-fold synonymous sites in the clade founder that are nucleotide $nt$, and $r_{nt_1\rightarrow nt_2}$ is the non-normalized rate of mutations from nucleotide $nt_1$ to $nt_2$ at 4-fold degenerate synonymous sites.

The expected number of mutations at each site (under neutrality) from the parental identity of $nt_1$ to some other identity of $nt_2$ is then simply $T \times r_{nt_1\rightarrow nt_2}$, which we will call the normalized rate for that clade.

We compute these expected numbers of mutations versus the actual numbers of mutations at each site, only considering actual mutations that are single nucleotide changes from the clade founder codon.

The expected and actual number of nucleotide mutation counts at each site are in [results/expected_vs_actual_mut_counts/expected_vs_actual_mut_counts.csv](results/expected_vs_actual_mut_counts/expected_vs_actual_mut_counts.csv).

### Computation of amino-acid mutation fitness effects
We then collapse the expected and actual counts for each amino-acid mutation, excluding the small number of sites that are in overlapping reading frames and are specified to exclude in [config.yaml](config.yaml) under the `gene_overlaps` key.
Note that in this aggregation, we exclude any amino acids with a codon for which at least one constituent nucleotide site is masked in `UShER`.

We estimate the fitness effect $\Delta f$ of each mutation as
$$ \Delta f = \log \left(\frac{n_{actual} + P}{n_{expected} + P}\right)$$
where $P$ is a pseudocount specified in [config.yaml](config.yaml) as `fitness_pseudocount`, and we are using the natural log.
So mutations with more counts than expected will have positive fitness effects, and those with less counts than expected will have negative fitness effects.

Note that these fitness effects will only be accurate of the number of expected counts is reasonably high.

The resulting fitness effect estimates are written to the following files:

 - [results/aa_fitness/aamut_fitness_all.csv](results/aa_fitness/aamut_fitness_all.csv): estimates aggregating across all clades.
 - [results/aa_fitness/aamut_fitness_by_clade.csv](results/aa_fitness/aamut_fitness_by_clade.csv): estimates for each individual clade.
 - [results/aa_fitness/aamut_fitness_by_subset.csv](results/aa_fitness/aamut_fitness_by_subset.csv): estimates splitting out by sequence subset (country).

Note also that the above files contain mutations in both numbering of the ORF1ab polypeptide and the nsp proteins contained within it.
The nsp protein mutations are a subset of the ORF1ab mutations, so if you examine both you would be double counting mutations.
The [config.yaml](config.yaml) file specifies the conversion from ORF1ab to nsp numbering.

Summaries are also plotted and are available

### Computation of amino-acid fitnesses
For each clade have estimated the change in fitness $\Delta f_{xy}$ caused by mutating a site from amino-acid $x$ to $y$, where $x$ is the amino acid in the clade founder sequence.
For each such mutation, we also have $n_{xy}$ which is the number of **expected** mutations from the clade founder amino acid $x$ to $y$.
These $n_{xy}$ values are important because they give some estimate of our "confidence" in the $\Delta f_{xy}$ values: if a mutation has high expected counts (large $n_{xy}$) then we can estimate the change in fitness caused by the mutation more accurately, and if $n_{xy}$ is small then the estimate will be much noisier.

However, we would like to aggregate the data across multiple clades to estimate amino-acid fitness values at a site under the assumption that these are constant across clades.
Now things get more complicated.
For instance, let's say at our site of interest, the clade founder amino acid is $x$ in one clade and $z$ in another clade.
For each clade we then have a set of $\Delta f_{xy}$ and $n_{xy}$ values for the first clade (where $y$ ranges over the 20 amino acids, including stop codon, that aren't $x$), and another set of up to 20 $\Delta f_{zy}$ and $n_{zy}$ values for the second clade (where $y$ ranges over the 20 amino acids that aren't $z$).

From these sets of mutation fitness changes, we'd like to estimate the fitness $f_x$ of each amino acid $x$, where the $f_x$ values satisfy $\Delta f_{xy} = f_y - f_x$ (in other words, a higher $f_x$ means higher fitness of that amino acid).
When there are multiple clades with different founder amino acids at the site, there is no guarantee that we can find $f_x$ values that precisely satisfy the above equation since there are more $\Delta f_{xy}$ values than $f_x$ values and the $\Delta f_{xy}$ values may have noise (and is some cases even real shifts among clades due to epistasis).
Nonetheless, we can try to find the $f_x$ values that come closest to satisfying the above equation.

First, we choose one amino acid to have a fitness value of zero, since the scale of the $f_x$ values is arbitrary and there are really only 20 unique parameters among the 21 $f_x$ values (there are 21 amino acids since we consider stops, but we only measure differences among them, not absolute values).
Typically if there was just one clade, we would set the wildtype value of $f_x = 0$ and then for mutations to all other amino acids $y$ we would simply have $f_y = \Delta f_{xy}$.
However, when there are multple clades with different founder amino acids, there is no longer a well defined "wildtype".
So we choose the most common **non-stop** parental amino-acid for the observed mutations and set that to zero.
In other words, we find $x$ that maximizes $\sum_y n_{xy}$ and set that $f_x$ value to zero.

Next, we choose the $f_x$ values that most closely match the measured mutation effects, weighting more strongly mutation effects with higher expected counts (since these should be more accurate).
Specifically, we define a loss function as
$$
L = \sum_x \sum_{y \ne x} n_{xy} \left(\Delta f_{xy} - \left[f_y - f_x\right]\right)^2
$$
where we ignore effects of synonymous mutations (the $x \ne y$ term in second summand) because we are only examining protein-level effects.
We then use numerical optimization to find the $f_x$ values that minimize that loss $L$.

Finally, we would still like to report an equivalent of the $n_{xy}$ values for the $\Delta f_{xy}$ values that give us some sense of how accurately we have estimated the fitness $f_x$ of each amino acid.
To do that, we tabulate $N_x = \sum_y \left(n_{xy} + n_{yx} \right)$ as the total number of mutations either from or to amino-acid $x$ as the "count" for the amino acid.
Amino acids with larger values of $N_x$ should have more accurate estimates of $f_x$.

The resulting amino-acid fitness values (aggregated across all clades) are in the following file:

 - [results/aa_fitness/aa_fitness.csv](results/aa_fitness/aa_fitness.csv)

The are also plotted in the heat maps at [https://jbloomlab.github.io/SARS2-mut-fitness](https://jbloomlab.github.io/SARS2-mut-fitness)

### Comparison to deep mutational scanning mutation effects
We compare the estimated fitness values to those extracted from a set of deep mutational scanning studies as specified under `dms_datasets` in [config.yaml](config.yaml).
The processed deep mutational scanning mutation effects are in `processed.csv` files in subdirectories of [./results/dms/](results/dms).

The correlation of the fitness estimates to the deep mutational scanning are plotted at [https://jbloomlab.github.io/SARS2-mut-fitness](https://jbloomlab.github.io/SARS2-mut-fitness)

### Non-terminal versus terminal counts
We compare the fitness effects of mutations to how often the mutation is observed on non-terminal versus terminal branches of the tree, and the mean log tip descendants sharing the mutations for each mutated branch.
See the plot linked at [https://jbloomlab.github.io/SARS2-mut-fitness](https://jbloomlab.github.io/SARS2-mut-fitness)

### Caveats of analysis
None of these are expected to seriously affect the accuracy of the current analysis, but they could become problematic if the same analysis is applied to substantially more diverged clades:

 - For computing the rates of mutations from 4-fold degenerate synonymous sites, we do not exclude sites / mutations to exclude from the computation of the sequence composition that is used to normalize the counts of different mutations to rates. This would become a problem if you start to specify a very large number of sites to exclude in [config.yaml](config.yaml), or if the sequences become highly diverged from the reference (as we exclude reversions to reference and their complement).
 
 - Multiple mutations in same codon on a branch are excluded from analysis. This is not expected to have much effect, as this is rare.

 - Indels are ignored.

 - The way that codon positions are assigned to identify synonymous sites will fail if there are non-in-frame indels.

 - Four-fold synonymous sites are identified in the clade founder, which could lead to mis-identification if seuqences in a clade become highly diverged from the founder.

 - Multiple mutations in the same codon in a clade can violate the assumptions about how sites are defined as synonymous, etc. For this reason they are excluded, and so we only include mutations that are from the clade founder amino acid identity.
 
 - We don't consider non-uniformity in mutation rate across the primary sequence.
 
 - If there are partial sequences in the tree (such that some sites are observed more than others), that is not accounted for.
