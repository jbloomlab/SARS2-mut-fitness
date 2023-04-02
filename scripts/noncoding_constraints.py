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
    window = np.ones(ws)
    return np.convolve(pos*ind, window, mode)/np.convolve(ind, window, mode), \
           np.convolve(values*ind, window, mode)/np.convolve(ind, window, mode)

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


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fitness", type=str, help="fitness file")
    parser.add_argument("--output", type=str, help="output file")
    args = parser.parse_args()

    fitness = pd.read_csv(args.fitness, sep=',')

    expected_count_cutoff = 20

    convolve_func = convolve_by_valid_entry
    windows = [3001, 11, 21, 51, 101]

    regions = {}
    peak = -0.7
    cutoff = -0.2
    region_intervals = np.zeros(29903)
    data = {}
    rois = {}
    fig, axs = plt.subplots(2,1,figsize=(18,6))
    ws0 = windows[0]
    cfit_by_pos = fitness.loc[(fitness['synonymous']|fitness['noncoding'])
                               &(fitness['expected_count']>expected_count_cutoff)
                            ].groupby('nt_site').mean().sort_index()


    values = cfit_by_pos['delta_fitness']
    data['position'] = np.arange(29903)
    data['raw'] = np.nan*np.ones_like(data['position'])
    positions = cfit_by_pos.index
    data['raw'][positions] = values #[values[p] if p in values else np.nan for p in data['position']]
    x0,y0 = convolve_func(positions, values,  ws0, mode='same')
    clade = 'all'
    for ws in windows[1:]:
        # ws_o_2 = int((ws0-ws)/2)
        ws_o_2 = int(ws/2)
        x,y = convolve_func(positions, values, ws)
        interpolator = interp1d(x,y,kind='linear', bounds_error=False)

        data[f'smooth_{ws}'] = [interpolator(p) for p in data['position']]
        axs[0].plot(x,y - y0[ws_o_2:-ws_o_2])
        axs[1].plot(x,y - y0[ws_o_2:-ws_o_2])

    axs[0].plot([0,30000], [0,0], lw=3, alpha=0.3, c='k')
    axs[1].plot([0,30000], [0,0], lw=3, alpha=0.3, c='k')
    add_genes(axs[0], 0.5, 0.8)
    add_genes(axs[1], 0.5, 0.8)
    axs[1].set_xlim(25000, 29700)
    # plt.savefig(f"results/figures/conservation_profile_{clade}.pdf")
    # d = pl.DataFrame(data)
    # d.write_csv(f"results/conservation_profile_{clade}.tsv", sep='\t')

