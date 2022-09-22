"""Get all sites in the GTF that correspond to a coding sequence."""


import pandas as pd


sites = (
    pd.read_csv(snakemake.input.gtf, sep="\t", header=None)
    .rename(columns={2: "feature_type", 3: "start", 4: "end"})
    [["feature_type", "start", "end"]]
    .query("feature_type == 'CDS'")
    .assign(
        site=lambda x: x.apply(lambda r: list(range(r["start"], r["end"] + 1)), axis=1),
    )
    .explode("site")
    ["site"]
    .drop_duplicates()
    .sort_values()
)

sites.to_csv(snakemake.output.csv, index=False)
