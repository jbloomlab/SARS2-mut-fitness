"""Save the amino acid fitness data in a JSON format along with metadata"""

import pandas as pd
import json

# Metadata for this dataset
citation = snakemake.params.citation
authors = snakemake.params.authors
source = snakemake.params.source
description = snakemake.params.description
tree = snakemake.wildcards.mat

# Load in the fitness data
aa_fitness = pd.read_csv(snakemake.input.fitness)

# Drop and rename columns
aa_fitness.drop(
    columns=["aa_differs_among_clade_founders", "subset_of_ORF1ab"], inplace=True
)
aa_fitness.rename(columns={"aa_site": "site", "aa": "mutant"}, inplace=True)

# Filter based on a minumum expected_count
min_expected_count = snakemake.params.min_expected_count
filtered_aa_fitness = aa_fitness.query("expected_count >= @min_expected_count")

# Group by gene and site and calculate the mean, min, max, median, and sum
filtered_aa_fitness_summary = (
    filtered_aa_fitness.query('mutant != "*"')  # Ignore stop codons in this calculation
    .groupby(["gene", "site"])["fitness"]
    .agg(["mean", "min", "max", "median", "sum"])
    .reset_index()
)

# Make a dictionary to dump into JSON
json_formatted_dict = {
    "citation": citation,
    "authors": authors,
    "source": source,
    "tree": tree,
    "description": description,
    "data": filtered_aa_fitness.to_dict("records"),
    "summary": filtered_aa_fitness_summary.to_dict("records"),
}

# Dump the dictionary into a JSON file
with open(snakemake.output.json, "w") as f:
    json.dump(json_formatted_dict, f, indent=4)
