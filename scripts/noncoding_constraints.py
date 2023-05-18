import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import argparse
from scipy.interpolate import interp1d
from known_conserved_noncoding import get_conserved_regions

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

def add_genes(ax, vmin, vmax, both=False, additional_regions = None, xscale=30000, region=None):
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
            if region is None or (f_start>=region[0] and f_start<region[1]):
                ax.text(f_start, vmin + offset, f.qualifiers['gene'][0])
                if both==True:
                    ax.text(f_stop, vmin + offset, f.qualifiers['gene'][0], ha='right')
            fi += 1
    if additional_regions is not None:
        for label, (f_start, f_stop) in sorted(additional_regions.items(), key=lambda x:x[1][0]):
            offset = (vmax - vmin) * (2 + fi%2)
            rect = patches.Rectangle((f_start, vmin+offset), f_stop - f_start, vmax-vmin, linewidth=1, edgecolor=f'C3', facecolor=f'C3')
            ax.add_patch(rect)
            if region is None or (f_start>=region[0] and f_start<region[1]):
                ax.text(f_stop+xscale/300, vmin + offset, label)
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
    parser.add_argument("--min_expected_count", type=int, help="min expected count")
    args = parser.parse_args()

    # read in tabular fitness file
    fitness = pd.read_csv(args.fitness, sep=',')
    # produce copies where ORFs that tolerate stop codons are marked as non-coding
    # give that aa substitutions don't matter there, we can use all sites to look for
    # non-coding constraints
    # don't include E and M regions
    known_conserved_positions, regions = get_conserved_regions(only_known=True)
    fitness['conserved_noncoding_known'] = fitness['nt_site'].apply(lambda x:bool(known_conserved_positions[x-1]))

    # include E and M regions
    conserved_noncoding_pos, regions = get_conserved_regions()
    fitness['conserved_noncoding'] = fitness["nt_site"].apply(lambda x: conserved_noncoding_pos[x-1])


    fitness_noORFs = erase_genes(fitness, ['ORF6', 'ORF7a', 'ORF7b','ORF8', 'ORF9b', 'ORF10'])

    expected_count_cutoff = args.min_expected_count
    fitness_measure = "delta_fitness"
    fs = 16

    # running window smoothing -- there is choice whether the smoothing window should run over a
    # range of ws adjacent nucleotides ('by_position') or over ws valid values without coding constraints
    convolve_mode = 'by_position'
    ws0 = 3000
    windows = [25,50]


    # subset table to position without plausible coding constraint, group by position, average
    cfit_by_pos = fitness_noORFs.loc[(fitness_noORFs['synonymous']|fitness_noORFs['noncoding'])
                               &(fitness_noORFs['expected_count']>expected_count_cutoff)
                            ].groupby('nt_site').mean().sort_index()

    values = cfit_by_pos[fitness_measure]
    positions = cfit_by_pos.index

    # set up figure
    fig = plt.figure(layout="constrained", figsize=(14,10))
    gs = GridSpec(4, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, :])
    ax3 = fig.add_subplot(gs[2:, :1])
    ax4 = fig.add_subplot(gs[2:, 1:])

    ## Make non-coding running average figure
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
    for ws in windows:
        ws_o_2 = int(ws/2)
        if convolve_mode=='by_position':
            x,y = convolve_by_position(data['position'], data['raw'], data['masked'], ws,  mode='valid')
        else:
            x,y = convolve_by_valid_entry(positions, values, ws)
        interpolator = interp1d(x,y,kind='linear', bounds_error=False)

        data[f'smooth_{ws}'] = [interpolator(p) for p in data['position']]
        ax1.plot(x,y, label=f"window {ws}nt") # - y0[ws_o_2:-ws_o_2])
        ax2.plot(x,y) # - y0[ws_o_2:-ws_o_2])

    ax1.legend(loc=3)
    ax1.plot([0,30000], [0,0], lw=3, alpha=0.3, c='k')
    ax2.plot([0,30000], [0,0], lw=3, alpha=0.3, c='k')
    add_genes(ax1, 0.5, 0.8, additional_regions=regions, xscale=30000, region=[0,25000])
    add_genes(ax2, 0.5, 0.8, additional_regions=regions, xscale=5000, region=[25000,30000])
    ax1.set_xlim(0, 29800)
    ax2.set_xlim(25000, 29800)
    ax1.set_ylabel('fitness')
    ax2.set_ylabel('fitness')
    ax2.set_xlabel('genome coordinate')
    ax1.text(-0.05,0.9, 'A', fontsize=fs*1.5, transform=ax1.transAxes)
    ax2.text(-0.05,0.9, 'B', fontsize=fs*1.5, transform=ax2.transAxes)

    ## Set up synonymous distribution figure
    all_ffold = fitness.loc[fitness['four_fold_degenerate']&(fitness['expected_count']>expected_count_cutoff)]
    conserved_ffold = fitness.loc[fitness['four_fold_degenerate']&fitness['conserved_noncoding']&(fitness['expected_count']>expected_count_cutoff)]
    other_ffold = fitness.loc[fitness['four_fold_degenerate']&(~fitness['conserved_noncoding'])&(fitness['expected_count']>expected_count_cutoff)]
    conserved_ffold_known = fitness.loc[fitness['four_fold_degenerate']&fitness['conserved_noncoding_known']&(fitness['expected_count']>expected_count_cutoff)]
    other_ffold_known = fitness.loc[fitness['four_fold_degenerate']&(~fitness['conserved_noncoding_known'])&(fitness['expected_count']>expected_count_cutoff)]

    ax3.hist(conserved_ffold_known[fitness_measure], bins = np.linspace(-6,2,20),
             label=f'known regions, n={len(conserved_ffold_known)}', alpha=0.5, density=True)
    ax3.hist(other_ffold_known[fitness_measure], bins = np.linspace(-6,2,20),
             label=f'other regions, n={len(other_ffold_known)}', alpha=0.5, density=True)
    ax3.legend()
    ax3.set_xlabel('fitness')
    ax3.text(-0.1,0.9, 'C', fontsize=fs*1.5, transform=ax3.transAxes)

    ax4.plot(sorted(all_ffold[fitness_measure]), np.arange(1,len(all_ffold)+1), label='all', c='C0')

    # ax4.plot(sorted(conserved_ffold[fitness_measure]), np.arange(1,len(conserved_ffold)+1),  label=f'known regions, n={len(conserved_ffold)}', c='C1', ls='--')
    # ax4.plot(sorted(other_ffold[fitness_measure]), np.arange(1,len(other_ffold)+1), label=f'other regions, n={len(other_ffold)}', c='C2', ls='--')

    ax4.plot(sorted(conserved_ffold_known[fitness_measure]), np.arange(1,len(conserved_ffold_known)+1),  label=f'known regions, n={len(conserved_ffold_known)}', c='C1')
    ax4.plot(sorted(other_ffold_known[fitness_measure]), np.arange(1,len(other_ffold_known)+1), label=f'other regions, n={len(other_ffold_known)}', c='C2')
    ax4.set_xlabel('fitness')
    ax4.set_ylabel('cumulative counts')
    ax4.set_yscale('log')
    ax4.legend()
    ax4.text(-0.1,0.9, 'D', fontsize=fs*1.5, transform=ax4.transAxes)

    plt.savefig(args.output, metadata={"creationDate": None})

