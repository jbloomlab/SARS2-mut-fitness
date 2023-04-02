from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def get_conserved_regions():
    TRS_motif = "ACGAAC"
    margin = 0
    ref = SeqIO.read('data/reference.gb', 'genbank')

    motifs = []
    pos=1
    while pos>0:
        pos = ref.seq.find(TRS_motif, pos)
        if pos>1:
            motifs.append(pos)
            pos += 1


    conserved_regions = []
    for mpos in motifs:
        conserved_regions.append((mpos-margin,mpos+len(TRS_motif)+margin))

    ## attenuator hairpin
    conserved_regions.append((ref.seq.find("ATGCTTCA"), ref.seq.find("CGTTTTT")))
    ## 3 stem pseudoknot
    conserved_regions.append((ref.seq.find("GCGGTGT"), ref.seq.find("TTTTGA", 13474)))
    ## center of E
    conserved_regions.append((26330,26360))

    conserved_vector = np.zeros(len(ref.seq), dtype=bool)
    for r in conserved_regions:
        conserved_vector[r[0]:r[1]] = True

    return conserved_vector


if __name__=="__main__":
    import pandas as pd

    fname = 'results/nt_fitness/ntmut_fitness_all.csv'
    min_expected_count = 20
    fitness = pd.read_csv(fname, sep=',')
    known_conserved_regions = get_conserved_regions()
    fitness['conserved_syn'] = fitness['nt_site'].apply(lambda x:bool(known_conserved_regions[x-1]))

    all_ffold = fitness.loc[fitness['four_fold_degenerate']&(fitness['expected_count']>min_expected_count)]
    conserved_ffold = fitness.loc[fitness['four_fold_degenerate']&fitness['conserved_syn']&(fitness['expected_count']>min_expected_count)]
    other_ffold = fitness.loc[fitness['four_fold_degenerate']&(~fitness['conserved_syn'])&(fitness['expected_count']>min_expected_count)]
    fitness_measure = 'delta_fitness'

    ## cumulative distribution of fitness effect of four-fold synonymous sites at known conserved regions vs other
    plt.figure()
    plt.plot(sorted(conserved_ffold[fitness_measure]), np.linspace(0,1,len(conserved_ffold)), label=f'known regions, n={len(conserved_ffold)}')
    plt.plot(sorted(other_ffold[fitness_measure]), np.linspace(0,1,len(other_ffold)), label=f'other regions, n={len(other_ffold)}')
    plt.legend()

    ## cumulative count of fitness effect of four-fold synonymous sites at known conserved regions vs other
    plt.figure()
    plt.plot(sorted(all_ffold[fitness_measure]), np.arange(0,len(all_ffold)), label='all')
    plt.plot(sorted(conserved_ffold[fitness_measure]), np.arange(0,len(conserved_ffold)),  label=f'known regions, n={len(conserved_ffold)}')
    plt.plot(sorted(other_ffold[fitness_measure]), np.arange(0,len(other_ffold)), label=f'other regions, n={len(other_ffold)}')
    plt.xlabel('fitness')
    plt.yscale('log')
    plt.legend()

    ## histogram of fitness effects
    plt.figure()
    plt.hist(all_ffold[fitness_measure], bins = np.linspace(-6,2,61), label='all', alpha=0.5)
    plt.hist(conserved_ffold[fitness_measure], bins = np.linspace(-6,2,61),  label=f'known regions, n={len(conserved_ffold)}', alpha=0.5)
    plt.hist(other_ffold[fitness_measure], bins = np.linspace(-6,2,61), label=f'other regions, n={len(other_ffold)}', alpha=0.5)
    plt.legend()
    plt.yscale('log')

    plt.figure()
    plt.hist(all_ffold[fitness_measure], bins = np.linspace(-6,2,61), label='all', alpha=0.5)
    plt.hist(conserved_ffold[fitness_measure], bins = np.linspace(-6,2,61),  label=f'known regions, n={len(conserved_ffold)}', alpha=0.5)
    plt.hist(other_ffold[fitness_measure], bins = np.linspace(-6,2,61), label=f'other regions, n={len(other_ffold)}', alpha=0.5)
    plt.legend()

    plt.figure()
    plt.hist(conserved_ffold[fitness_measure], bins = np.linspace(-6,2,16), label=f'known regions (n={len(conserved_ffold)})', alpha=0.5, density=True)
    plt.hist(other_ffold[fitness_measure], bins = np.linspace(-6,2,16), label=f'other regions, n={len(other_ffold)}', alpha=0.5, density=True)
    plt.xlabel('fitness')
    plt.ylabel('density')
    plt.legend()
