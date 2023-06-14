"""Save the amino acid fitness data in a JSON format along with metadata"""

import pandas as pd
import json
import gzip


# Metadata for this dataset
citation = snakemake.params.citation
authors = snakemake.params.authors
source = snakemake.params.source
description = snakemake.params.description
tree = snakemake.wildcards.mat

# Load in the fitness data
aa_fitness = pd.read_csv(snakemake.input.aa_fitness)

# Drop and rename columns
aa_fitness.drop(
    columns=["aa_differs_among_clade_founders", "subset_of_ORF1ab"], inplace=True
)
aa_fitness.rename(columns={"aa_site": "site", "aa": "mutant"}, inplace=True)

# Join the wildtype amino acids from Nextstrain clade 19A into the dataframe
clade_founder_aas = pd.read_csv(snakemake.input.clade_founder_aas)
clade_founder_aas_19A = (
    clade_founder_aas.query("clade == '19A (B)'")
    .drop(columns=["clade"])
    .rename(columns={"amino acid": "wildtype"})
)
aa_fitness = aa_fitness.merge(
    clade_founder_aas_19A, how="left", on=["gene", "site"], validate="many_to_one"
)

# Filter based on a minumum expected_count
min_expected_count = snakemake.params.min_expected_count
filtered_aa_fitness = aa_fitness.query("expected_count >= @min_expected_count")

# Group by gene and site and calculate the mean, min, max, median, and sum
summary_stats = ["mean", "min", "max", "median", "sum"]
filtered_aa_fitness_summary = (
    filtered_aa_fitness.query('mutant != "*"')  # Ignore stop codons in this calculation
    .groupby(["gene", "site", "wildtype"])["fitness"]
    .agg(summary_stats)
    .reset_index()
)
filtered_aa_fitness_summary.columns = [
    col + "_fitness" if col in summary_stats else col
    for col in filtered_aa_fitness_summary.columns
]

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
with open(snakemake.output.aa_fitness_json, "w") as f:
    json.dump(json_formatted_dict, f, indent=4)

# Compress the JSON file
with open(snakemake.output.aa_fitness_json, "rb") as f_in:
    with gzip.open(snakemake.output.aa_fitness_json_gz, "wb") as f_out:
        f_out.writelines(f_in)
