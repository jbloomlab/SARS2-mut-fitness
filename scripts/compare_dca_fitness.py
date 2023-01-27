import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

remap_genes = {"Spike":"S", "ORF1a":"ORF1ab", "ORF1b":"ORF1ab",
                "Nucleocapsid":"N", "Membrane":"M", "Envelope":"E"}

def load_dca():
    dca = pd.read_csv(snakemake.input.dca)
    dca_mutability = {}
    for i, row in dca.iterrows():
        dca_mutability[(remap_genes.get(row["protein"], row["protein"]), row["position_protein"])] = row["mutability_score(DCA)"]
    return dca_mutability



def load_fitness():
    return pd.read_csv(snakemake.input.fitness)


if __name__ == "__main__":
    dca_mutability = load_dca()
    fitness = load_fitness()
    fitness = fitness.loc[fitness.expected_count>30]

    average_fitness = fitness.groupby(["gene", "aa_site"])['delta_fitness'].mean()

    linked_dca_fitness = defaultdict(dict)
    for i, row in average_fitness.iteritems():
        gene, aa_site = i
        if (gene, aa_site) in dca_mutability and not pd.isna(dca_mutability[(gene, aa_site)]):
            linked_dca_fitness[gene][aa_site] = (dca_mutability.get((gene, aa_site), np.nan), row)

    fig, axes = plt.subplots(2, 5, figsize=(12, 5))
    for i, gene in enumerate(linked_dca_fitness):
        ax = axes.ravel()[i]
        d = np.array(list(linked_dca_fitness[gene].values()))
        corr = np.corrcoef(d.T)[0,1]
        ax.scatter(d[:,0], d[:,1], s=10, alpha=0.3)
        ax.set_title(f"{gene} (r = {corr:.2f})")
        ax.set_xlabel('DCA mutability')
        ax.set_ylabel('mean fitness effect')
    fig.tight_layout()
    
    fig.savefig(snakemake.output.plot)
