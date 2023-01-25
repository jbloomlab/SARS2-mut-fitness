import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

remap_genes = {"Spike":"S", "ORF1a":"ORF1ab", "ORF1b":"ORF1ab",
                "Nucleocapsid":"N", "Membrane":"M", "Envelope":"E"}

def load_dca():
    dca = pd.read_csv("data/dca_mutability.csv")
    dca_mutability = {}
    for i, row in dca.iterrows():
        dca_mutability[(remap_genes.get(row["protein"], row["protein"]), row["position_protein"])] = row["mutability_score(DCA)"]
    return dca_mutability



def load_fitness():
    return pd.read_csv("results/aa_fitness/aamut_fitness_all.csv")


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


    for gene in linked_dca_fitness:
        d = np.array(list(linked_dca_fitness[gene].values()))
        print(gene, np.corrcoef(d.T)[0,1])
        plt.figure()
        plt.scatter(d[:,0], d[:,1])
        plt.title(gene)
        plt.xlabel('dca')
        plt.ylabel('fitness')
