import itertools

import pandas as pd

import yaml

nts = ["A", "C", "G", "T"]

sites_to_exclude = pd.DataFrame(
    [
        (clade, site, f"{nt1}{site}{nt2}", False)
        for nt1, site, nt2, clade in itertools.product(
            nts,
            (
                snakemake.params.sites_to_exclude
                + pd.read_csv(snakemake.input.site_mask)["site"].tolist()
            ),
            nts,
            snakemake.params.clades,
        )
    ],
    columns=["clade", "site", "mutation", "masked_in_usher"],
).drop_duplicates()

if snakemake.params.exclude_ref_to_founder_muts:
    muts_to_exclude = pd.concat(
        [
            pd.read_csv(f).assign(clade=clade)
            for f, clade in zip(
                snakemake.input.muts_to_exclude, snakemake.params.clades
            )
        ]
    ).assign(masked_in_usher=False)
    to_exclude = pd.concat([sites_to_exclude, muts_to_exclude])
else:
    to_exclude = sites_to_exclude

with open(snakemake.input.usher_masked_sites) as f:
    usher_masked_sites = yaml.safe_load(f)
for mask_dict in usher_masked_sites.values():
    mask_df = pd.DataFrame(
        [
            (clade, site, f"{nt1}{site}{nt2}", True)
            for nt1, site, nt2, clade in itertools.product(
                nts,
                mask_dict["sites"],
                nts,
                mask_dict["clades"],
            )
        ],
        columns=["clade", "site", "mutation", "masked_in_usher"],
    )
    to_exclude = pd.concat([to_exclude, mask_df])

(
    to_exclude
    .groupby(["clade", "site", "mutation"], as_index=False)
    .aggregate(masked_in_usher=pd.NamedAgg("masked_in_usher", "any"))
    .sort_values(["clade", "site", "mutation"])
    .to_csv(snakemake.output.csv, index=False)
)
