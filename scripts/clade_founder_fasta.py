"""Get clade founder FASTA from JSON, indels ignored."""


import json

import Bio.SeqIO

import pandas as pd


input_neher_json = snakemake.input.neher_json
input_roemer_json = snakemake.input.roemer_json
roemer_nextstrain_to_pango = snakemake.params.roemer_nextstrain_to_pango
ref_fasta = snakemake.input.ref_fasta
clade = snakemake.wildcards.clade


print(f"Getting founder fasta for {clade=}")

with open(input_neher_json) as f:
    neher_founder_muts = json.load(f)
with open(input_roemer_json) as f:
    roemer_founder_muts = json.load(f)

if clade in neher_founder_muts:
    muts = neher_founder_muts[clade]["nuc"]
    assert clade not in roemer_nextstrain_to_pango, f"{clade=} in Neher and Roemer"

elif clade in roemer_nextstrain_to_pango:
    pango = roemer_nextstrain_to_pango[clade]
    muts = roemer_founder_muts[pango]["nucSubstitutions"]
    assert roemer_founder_muts[pango]["nextstrainClade"] == clade

else:
    raise ValueError(f"no founder for {clade=}")

print(f"Clade has the following {len(muts)} nucleotide mutations:\n{muts}")

ref = list(str(Bio.SeqIO.read(ref_fasta, "fasta").seq))

ref_to_founder_muts = []
for m in muts:
    wt, site, mut = m[0], int(m[1:-1]), m[-1]
    assert ref[site - 1] == wt, "{m=}, {ref[site - 1]=}"
    ref[site - 1] = mut
    ref_to_founder_muts += [(site, m), (site, f"{mut}{site}{wt}")]

ref = "".join(ref)

with open(snakemake.output.fasta, "w") as f:
    f.write(f">{clade}\n{ref}\n")

pd.DataFrame(ref_to_founder_muts, columns=["site", "mutation"]).to_csv(
    snakemake.output.muts,
    index=False,
)
