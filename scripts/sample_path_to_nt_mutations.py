"""Convert `matUtils` sample paths to mutations leading to each node."""


import pandas as pd


# get nodes and mutation list for each node
nodes_and_muts = (
    pd.read_csv(snakemake.input.tsv, sep="\t", names=["tip", "path"])
    .assign(node_paths=lambda x: x["path"].str.split())
    .drop(columns=["tip", "path"])
    .explode("node_paths")
    .query("node_paths.notnull()")
    .assign(
        node_id=lambda x: x["node_paths"].str.split(":").str[0],
        nt_mutations=lambda x: x["node_paths"].str.split(":").str[1],
    )
    .drop(columns="node_paths")
)

# ensure each node has a unique muation path
assert nodes_and_muts["node_id"].nunique() == len(nodes_and_muts.drop_duplicates())

# get just unique entry per node and clean up
nodes_and_muts = (
    nodes_and_muts
    .drop_duplicates()
    .assign(nt_mutations=lambda x: x["nt_mutations"].str.replace(",", ";"))
)

nodes_and_muts.to_csv(snakemake.output.csv, index=False)
