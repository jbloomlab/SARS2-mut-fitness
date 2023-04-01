import os

import Bio.SeqIO

import pandas as pd


# codons that start with these two nucleotides are 4-fold degenerate
codons_4fold = {
    "TC",  # serine
    "CT",  # leucine
    "CC",  # proline
    "CG",  # arginine
    "AC",  # threonine
    "GT",  # valine
    "GC",  # alanine
    "GG",  # glycine
}

sites = (
    pd.read_csv(snakemake.input.coding_sites)
    .sort_values("site")
    .assign(
        codon_position=lambda x: x["codon_position"].str.split(";"),
        codon_site=lambda x: x["codon_site"].str.split(";"),
        gene=lambda x: x["gene"].str.split(";"),
    )
    .explode(["codon_position", "codon_site", "gene"])
)


def get_codon(row):
    if row["codon_position"] == "noncoding":
        return "noncoding"
    r = row["site"]
    codon_pos = int(row["codon_position"])
    return f"{seq[r - codon_pos]}{seq[r - codon_pos + 1]}{seq[r - codon_pos + 2]}"


dfs = []
for f_name in snakemake.input.fastas:
    seq = str(Bio.SeqIO.read(f_name, "fasta").seq)
    name = os.path.splitext(os.path.basename(f_name))[0]
    dfs.append(
        sites.assign(
            clade=name,
            nt=lambda x: x["site"].map(lambda r: seq[r - 1]),  # noqa: B023
            codon=lambda x: x.apply(get_codon, axis=1),
            four_fold_degenerate=lambda x: (
                (x["codon_position"] == "3") & x["codon"].str[:2].isin(codons_4fold)
            ),
        )
        .groupby(["clade", "site", "nt"], as_index=False)
        .aggregate(
            gene=pd.NamedAgg("gene", lambda s: ";".join(s)),
            codon=pd.NamedAgg("codon", lambda s: ";".join(s)),
            codon_position=pd.NamedAgg(
                "codon_position", lambda s: ";".join(map(str, s))
            ),
            codon_site=pd.NamedAgg("codon_site", lambda s: ";".join(s)),
            four_fold_degenerate=pd.NamedAgg("four_fold_degenerate", "all"),
        )
    )

pd.concat(dfs).to_csv(snakemake.output.csv, index=False)
