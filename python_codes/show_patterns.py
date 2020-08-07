import sys
import chromosight.utils.contacts_map as cuc
import pandas as pd
import cooler
import numpy as np
import matplotlib.pyplot as plt

cool = 'GSE63525_GM12878_insitu_primary.mcool::/resolutions/10000'
region = "2:190,800,000-191,900,000"
#region = "2:150,800,000-200,900,000"


c = cooler.Cooler(cool)
chromo = pd.read_csv('chromosight/chromosight_loops_small_GSE63525.tsv', sep='\t')
hiccups = pd.read_csv('hiccups/GSE63525_GM12878_primary_HiCCUPS_looplist.txt', sep='\t')
hicexplorer = pd.read_csv('hicexplorer/hicexplorer_loops.tsv', sep='\t', header=None, usecols=[3, 4, 0, 1], names=['chr1', 'x', 'chr2', 'y'])
cooltools = pd.read_csv('cooltools/cooltools_loops.tsv.postproc', sep='\t', usecols=[0, 1, 3, 4])

hg = cuc.HicGenome(cool)

hiccups['bin1'] = hg.coords_to_bins(hiccups.loc[:, ['chr1', 'x1']].rename(columns={'chr1': 'chrom', 'x1': 'pos'}))
hiccups['bin2'] = hg.coords_to_bins(hiccups.loc[:, ['chr2', 'y1']].rename(columns={'chr2': 'chrom', 'y1': 'pos'}))

hicexplorer['bin1'] = hg.coords_to_bins(hicexplorer.loc[:, ['chr1', 'x']].rename(columns={'chr1': 'chrom', 'x': 'pos'}))
hicexplorer['bin2'] = hg.coords_to_bins(hicexplorer.loc[:, ['chr2', 'y']].rename(columns={'chr2': 'chrom', 'y': 'pos'}))

cooltools['bin1'] = hg.coords_to_bins(cooltools.loc[:, ['chrom1', 'start1']].rename(columns={'chrom1': 'chrom', 'start1': 'pos'}))
cooltools['bin2'] = hg.coords_to_bins(cooltools.loc[:, ['chrom2', 'start2']].rename(columns={'chrom2': 'chrom', 'start2': 'pos'}))

mat = c.matrix().fetch(region)
mat = mat**0.2
mat[mat == 0] = np.nan

s1, e1 = c.extent(region)
s2, e2 = s1, e1

fig, ax = plt.subplots(1, 4, sharex=True, sharey=True)
plt.suptitle('STAT1, chr2:191,829,084-191,885,686')

for i, pat, label in zip(
    range(4),
    [chromo, hiccups, hicexplorer, cooltools],
    ['chromosight', 'hiccups', 'hicexplorer', 'cooltools'],
):
    pat = pat.loc[
        (pat.bin1 > s1) & (pat.bin1 < e1) & (pat.bin2 > s2) & (pat.bin2 < e2),
        :,
    ]
    ax[i].imshow(
        mat,
        cmap="afmhot_r",
        vmax=np.percentile(mat[~np.isnan(mat)], 99),
    )
    ax[i].scatter(
        pat.bin1 - s1,
        pat.bin2 - s2,
        facecolors="none",
        edgecolors='blue',
        marker='o',
        label=label
    )
    ax[i].set_title(label)

plt.show()
