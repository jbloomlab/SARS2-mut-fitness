"""Get all sites in the GTF that correspond to a coding sequence."""


import Bio.SeqIO

import pandas as pd


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
        gene=lambda x: x["attribute"].str.extract(r'^gene_id "(\w+)"'),
        # handle shifts for multiple CDSs for ORF1ab
        genelength=lambda x: x["end"] - x["start"] + 1,
        shift=lambda x: x.groupby("gene")["genelength"].shift(fill_value=0),
    )
    .drop(columns=["seqName", "source", "score", "feature", "attribute"])
)

assert all(cds_features["frame"] == 0), "non frame 0 CDSs\n" + str(cds_features)
assert all(cds_features["strand"] == "+"), "non + sense CDS\n" + str(cds_features)

# get sites with codon position assigned
sites = (
    cds_features.assign(
        site=lambda x: x.apply(lambda r: list(range(r["start"], r["end"] + 1)), axis=1),
    )
    .explode("site", ignore_index=True)
    .assign(
        codon_position=lambda x: x.apply(
            lambda r: (r["site"] - r["start"]) % 3 + 1,
            axis=1,
        ),
        codon_site=lambda x: x.apply(
            lambda r: (r["site"] - r["start"] + r["shift"]) // 3 + 1,
            axis=1,
        ),
    )
    .drop(columns=["strand", "frame", "start", "end"])
    .groupby("site", as_index=False)
    .aggregate(
        codon_position=pd.NamedAgg(
            "codon_position",
            lambda s: ";".join(map(str, s.values)),
        ),
        codon_site=pd.NamedAgg("codon_site", lambda s: ";".join(map(str, s.values))),
        gene=pd.NamedAgg("gene", lambda s: ";".join(s)),
    )
)

# add any non-coding sites
length = len(Bio.SeqIO.read(snakemake.input.fasta, "fasta"))
noncoding_sites = sorted(set(range(1, length + 1)) - set(sites["site"]))

sites = pd.concat([
    sites,
    pd.DataFrame(
        {
            "site": noncoding_sites,
            "codon_position": "noncoding",
            "codon_site": "noncoding",
            "gene": "noncoding",
        }
    ),
]).sort_values("site")

sites.to_csv(snakemake.output.csv, index=False)
