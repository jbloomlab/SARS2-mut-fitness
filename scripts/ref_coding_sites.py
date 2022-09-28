"""Get all sites in the GTF that correspond to a coding sequence."""


import pandas as pd

import re


# get CDS features
cds_features = (
    pd.read_csv(
        snakemake.input.gtf,
        sep="\t",
        header=None,
        names=[
            "seqName",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )
    .query("feature == 'CDS'")
    .assign(
        frame=lambda x: x["frame"].astype(int),
        gene=lambda x: x["attribute"].str.extract('^gene_id "(\w+)"'),
    )
    .drop(columns=["seqName", "source", "score", "feature", "attribute"])
)

assert all(cds_features["frame"] == 0), "non frame 0 CDSs\n" + str(cds_features)
assert all(cds_features["strand"] == "+"), "non + sense CDS\n" + str(cds_features)

# get sites with codon position assigned
sites = (
    cds_features
    .assign(
        site=lambda x: x.apply(lambda r: list(range(r["start"], r["end"] + 1)), axis=1),
    )
    .explode("site", ignore_index=True)
    .assign(
        codon_position=lambda x: x.apply(
            lambda r: (r["site"] - r["start"]) % 3 + 1, axis=1,
        ),
    )
    .drop(columns=["strand", "frame", "start", "end"])
    .groupby("site", as_index=False)
    .aggregate(
        codon_position=pd.NamedAgg(
            "codon_position", lambda s: ";".join(map(str, s.values)),
        ),
        gene=pd.NamedAgg("gene", lambda s: ";".join(s)),
    )
)

sites.to_csv(snakemake.output.csv, index=False)
