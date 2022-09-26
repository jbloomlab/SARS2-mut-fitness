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

sites = pd.read_csv(snakemake.input.coding_sites).sort_values("site")

dfs = []
for f_name in [snakemake.input.ref_fasta, *snakemake.input.fastas]:
    seq = str(Bio.SeqIO.read(f_name, "fasta").seq)
    name = os.path.splitext(os.path.basename(f_name))[0]
    dfs.append(
        sites
        .assign(
            clade=name,
            nt=lambda x: x["site"].map(lambda r: seq[r - 1]),
            preceding_2nt=lambda x: x["site"].map(lambda r: seq[r - 3: r - 1]),
            four_fold_degenerate=lambda x: (
                (x["codon_position"].astype(str) == "3")
                & x["preceding_2nt"].isin(codons_4fold)
            )
        )
    )

pd.concat(dfs).to_csv(snakemake.output.csv, index=False)
