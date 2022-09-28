"""Get samples from mutation-annotated tree and annotate clade and date."""


import subprocess

import pandas as pd


subprocess.run(
    ["matUtils", "summary", "-i", snakemake.input.mat, "-C", snakemake.output.csv]
)

df = (
    pd.read_csv(snakemake.output.csv, sep="\t")
    .rename(columns={"annotation_1": "nextstrain_clade", "annotation_2": "pango_clade"})
    .assign(
        nextstrain_clade=lambda x: x["nextstrain_clade"].str.split().str[0],
        date=lambda x: x["sample"].str.split("|").str[-1],
    )
)

df.to_csv(snakemake.output.csv, index=False)

clade_counts = (
    df.groupby("nextstrain_clade", as_index=False)
    .aggregate(count=pd.NamedAgg("sample", "count"))
    .assign(
        adequate_sample_counts=lambda x: x["count"]
        >= snakemake.params.min_clade_samples
    )
)

clade_counts.to_csv(snakemake.output.clade_counts, index=False)
