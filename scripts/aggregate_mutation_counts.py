import os

import pandas as pd


counts = pd.concat(
    [
        pd.read_csv(f).assign(
            clade=os.path.splitext(os.path.basename(f))[0].split("_")[0],
            subset=os.path.splitext(os.path.basename(f))[0].split("_")[1],
        )
        for f in snakemake.input.counts
    ]
)

clade_founder_nts = pd.read_csv(snakemake.input.clade_founder_nts).rename(
    columns={"site": "nt_site", "gene": "protein", "nt": "clade_founder_nt"},
)

# check that each site corresponds to same clade founder and protein
# if they correspond, the following merge should be one-to-one with no NaN
merge_cols = ["nt_site", "protein", "clade", "clade_founder_nt"]
merge_check = (
    counts[merge_cols]
    .drop_duplicates()
    .merge(clade_founder_nts, validate="one_to_one", how="left")
)
assert merge_check.notnull().all().all(), merge_check.notnull().all()

counts_and_info = counts.merge(
    clade_founder_nts,
    on=merge_cols,
    how="left",
    validate="many_to_one",
)

counts_and_info.to_csv(snakemake.output.csv, index=False)
