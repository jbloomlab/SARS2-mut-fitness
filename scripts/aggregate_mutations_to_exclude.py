import itertools

import pandas as pd


sites_to_exclude = pd.DataFrame(
    [
        (clade, site, f"{nt1}{site}{nt2}")
        for nt1, site, nt2, clade in itertools.product(
            ["A", "C", "G", "T"],
            snakemake.params.sites_to_exclude,
            ["A", "C", "G", "T"],
            snakemake.params.clades,
        )
    ],
    columns=["clade", "site", "mutation"],
)

if snakemake.params.exclude_ref_to_founder_muts:
    muts_to_exclude = pd.concat(
        [
            pd.read_csv(f).assign(clade=clade)
            for f, clade in zip(
                snakemake.input.muts_to_exclude, snakemake.params.clades
            )
        ]
    )
    to_exclude = pd.concat([sites_to_exclude, muts_to_exclude])
else:
    to_exclude = sites_to_exclude


(
    to_exclude.sort_values(["clade", "site", "mutation"]).to_csv(
        snakemake.output.csv, index=False
    )
)
