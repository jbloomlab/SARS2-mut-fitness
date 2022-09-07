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
        date=lambda x: x["sample"].str.split("|").str[-1],
        valid_date=lambda x: x["date"].str.fullmatch("\d{4}\-\d{2}\-\d{2}"),
    )
)

df.to_csv(snakemake.output.csv, index=False)
