"""Get clade founder FASTA from JSON, indels ignored."""


import json

import Bio.SeqIO

import pandas as pd


input_json = snakemake.input.json
ref_fasta = snakemake.input.ref_fasta
clade = snakemake.wildcards.clade


print(f"Getting founder fasta for {clade=}")

with open(input_json) as f:
    founder_muts = json.load(f)

muts = founder_muts[clade]["nuc"]

print(f"Clade has the following {len(muts)} nucleotide mutations:\n{muts}")

ref = list(str(Bio.SeqIO.read(ref_fasta, "fasta").seq))

ref_to_founder_muts = []
for m in muts:
    wt, site, mut = m[0], int(m[1: -1]), m[-1]
    assert ref[site - 1] == wt, "{m=}, {ref[site - 1]=}"
    ref[site - 1] = mut
    ref_to_founder_muts += [(site, m), (site, f"{mut}{site}{wt}")]

ref = "".join(ref)

with open(snakemake.output.fasta, "w") as f:
    f.write(f">{clade}\n{ref}\n")

pd.DataFrame(ref_to_founder_muts, columns=["site", "mutation"]).to_csv(
    snakemake.output.muts, index=False,
)
