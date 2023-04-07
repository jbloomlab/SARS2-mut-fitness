import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.interpolate import interp1d

def convolve_by_valid_entry(position, values, ws, mode='valid'):
    window = np.ones(ws)/ws
    return np.convolve(position, window, mode), np.convolve(values, window, mode)

def convolve_by_position(position, values, ind, ws, mode='valid'):
    pos = position
    tmp_values = np.copy(values)
    tmp_values[~ind]=0
    window = np.ones(ws)
    return np.convolve(pos*ind, window, mode)/np.convolve(ind, window, mode), \
           np.convolve(tmp_values*ind, window, mode)/np.convolve(ind, window, mode)

def add_genes(ax, vmin, vmax, both=False):
    import matplotlib.patches as patches
    from Bio import SeqIO
    ref = SeqIO.read('data/reference.gb', 'genbank')
    fi = 0
    for f in ref.features:
        if f.type=='CDS':
            offset = (vmax - vmin)/3 * (fi%2)
            f_start = f.location.start.position
            f_stop = f.location.end.position
            rect = patches.Rectangle((f_start, vmin+offset), f_stop - f_start, vmax-vmin, linewidth=1, edgecolor=f'C{fi}', facecolor=f'C{fi}')
            ax.add_patch(rect)
            ax.text(f_start, vmin + offset, f.qualifiers['gene'][0])
            if both==True:
                ax.text(f_stop, vmin + offset, f.qualifiers['gene'][0], ha='right')
            fi += 1


def erase_genes(fitness, genes):
    will_be_noncoding = fitness['gene'].apply(lambda x: x in genes)
    fitness.loc[will_be_noncoding,"noncoding"] = True
    fitness.loc[will_be_noncoding,"synonymous"] = False
    fitness.loc[will_be_noncoding,"four_fold_degenerate"] = False
    if 'ORF9b' in genes:
        from Bio import SeqIO, Seq
        ref = SeqIO.read('data/reference.gb', 'genbank')
        will_be_syn = []
        will_be_ffold = []
        N_start = 28273
        for row in fitness.itertuples():
            if row.gene=='N;ORF9b':
                a, pos, d = row.nt_mutation[0], row.nt_site-1, row.nt_mutation[-1]
                frame = (pos-N_start)%3
                codon = ref.seq[pos-frame:pos-frame+3]
                mut_codon = codon[:frame] + d + codon[frame+1:]
                will_be_syn.append(Seq.translate(codon)==Seq.translate(mut_codon))
                if frame==2 and all([Seq.translate(codon)==Seq.translate(codon[:2]+nuc) for nuc in 'ACGT']):
                    will_be_ffold.append(True)
                else:
                    will_be_ffold.append(False)
            else:
                will_be_syn.append(False)
                will_be_ffold.append(False)

    fitness.loc[will_be_syn,"synonymous"] = True
    fitness.loc[will_be_ffold,"four_fold_degenerate"] = True

    return fitness

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fitness", type=str, help="fitness file")
    parser.add_argument("--output", type=str, help="output file")
    args = parser.parse_args()

    # read in tabular fitness file
    fitness = pd.read_csv(args.fitness, sep=',')
    # produce copies where ORFs that tolerate stop codons are marked as non-coding
    # give that aa substitutions don't matter there, we can use all sites to look for
    # non-coding constraints
    fitness_noORFs = erase_genes(fitness, ['ORF6', 'ORF7a', 'ORF7b','ORF8', 'ORF9b', 'ORF10'])

    expected_count_cutoff = 20

    # running window smoothing -- there is choice whether the smoothing window should run over a
    # range of ws adjacent nucleotides ('by_position') or over ws valid values without coding constraints
    convolve_mode = 'by_position'
    ws0 = 3000
    windows = [24,46]


    # subset table to position without plausible coding constraint, group by position, average
    cfit_by_pos = fitness_noORFs.loc[(fitness_noORFs['synonymous']|fitness_noORFs['noncoding'])
                               &(fitness_noORFs['expected_count']>expected_count_cutoff)
                            ].groupby('nt_site').mean().sort_index()

    values = cfit_by_pos['delta_fitness']
    positions = cfit_by_pos.index

    data = {}
    data['position'] = np.arange(29903)
    data['raw'] = np.nan*np.ones_like(data['position'])
    data['raw'][positions] = values
    data['masked'] = ~np.isnan(data['raw'])
    # generate large scale smoothed estimates as a baseline
    if convolve_mode=='by_position':
        x0,y0 = convolve_by_position(data['position'], data['raw'], data['masked'], ws0, mode='same')
    else:
        x0,y0 = convolve_by_valid_entry(positions, values,  ws0, mode='same')

    # loop over smoothing window sizes and plot the smoothed non-coding fitness estimates
    fig, axs = plt.subplots(2,1,figsize=(18,6))
    for ws in windows:
        ws_o_2 = int(ws/2)
        if convolve_mode=='by_position':
            x,y = convolve_by_position(data['position'], data['raw'], data['masked'], ws,  mode='valid')
        else:
            x,y = convolve_by_valid_entry(positions, values, ws)
        interpolator = interp1d(x,y,kind='linear', bounds_error=False)

        data[f'smooth_{ws}'] = [interpolator(p) for p in data['position']]
        axs[0].plot(x,y) # - y0[ws_o_2:-ws_o_2])
        axs[1].plot(x,y) # - y0[ws_o_2:-ws_o_2])

    axs[0].plot([0,30000], [0,0], lw=3, alpha=0.3, c='k')
    axs[1].plot([0,30000], [0,0], lw=3, alpha=0.3, c='k')
    add_genes(axs[0], 0.5, 0.8)
    add_genes(axs[1], 0.5, 0.8)
    axs[1].set_xlim(25000, 29700)
    plt.savefig("_noncoding_conservation.pdf")

