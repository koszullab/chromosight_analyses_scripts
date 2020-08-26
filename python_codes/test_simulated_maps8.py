# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:27:26 2018
7: reprise feb 2020
Notably we try other pattern of borders (the generic ones) to make propensity maps.
"""
import random
import itertools
import os
import sys
import time
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import statsmodels.api as sm
import cooler
import chromosight.utils.preprocessing as cup
import chromosight.kernels as ck

# Parse CLI arguments
# Path to the cool file with experimental data
clr = cooler.Cooler(sys.argv[1])
# Path to positions of borders in the experimental matrix (in bins)
borders_pos = np.loadtxt(sys.argv[2])


lowess = sm.nonparametric.lowess

# Experimental data:
m = clr.matrix(balance=True, sparse=False)[:]
bins = clr.bins[:]
bad_bins = bins.loc[np.isnan(bins.weight)].index.values
Nreads = sum(m)
Nrealisations = 10
n1 = shape(m)[0]  # size of the chromosome contact map

# Normalisation:


# Load patterns
# Generic pattern:
RATIO_CIDs = ck.borders["kernels"][0]
RATIO_CIDs = RATIO_CIDs ** 0.2  # to attenuate a bit
RATIO_LOOPS = ck.loops["kernels"][0]
RATIO_LOOPS = RATIO_LOOPS ** 0.7  # to attenuate a bit

# Remove borders from experimental data overlapping empty bins:
chrm = 5
border_pos = list(borders_pos)
for b in borders_pos.copy():
    # If any bad (empty) bin is closer to a border than the radius of the
    # pattern template, drop it
    if np.any(np.abs(bad_bins - b) < max(ratio_cids.shape)) in b:
        borders_pos.remove(b)
# To convert the 1D vetor into matrice object:
MAT_INDICES = np.zeros((shape(m)))
compt = 0
VECT_COMPT = []
for i, j in itertools.product(range(n1), range(n1)):
    MAT_INDICES[i, j] = compt
    VECT_COMPT.append((i, j))
    compt += 1

#  Genomic distance law from experiemental data:
s = shape(matscn1)[0]
prob_d = distance_law_human.dist_law(matscn1)
prob_d = prob_d / sum(prob_d)

# LOWESS Fit:
z = lowess(prob_d[:, 0], range(0, len(prob_d)), frac=1.0 / 15)
prob_d_lowess = z[:, 1]

# TADs size distribution:
list_size = []
for i in range(1, len(b)):
    bi = b[i] - b[i - 1]
    list_size.append(bi)

s = shape(matscn1)[0]

# Genomic distance law:
MAT_GENO = np.zeros((shape(m)))
for i in range(0, s):
    for j in range(0, s):
        MAT_GENO[i, j] = prob_d_lowess[abs(j - i)]

if len(MAT_GENO[MAT_GENO < 0]) > 0:
    raise ValueError(
        f"Presence of {len(mat_geno[mat_geno < 0])} negative elements in the "
        "distance law. Try changing Loess parameters."
    )
MAT_GENO[MAT_GENO < 0] = 0


# ------------------------------------------------------------------------------
#  RANDOM MAPs GENERATION
# ------------------------------------------------------------------------------
# NOTE: Send to output folder
t1 = time.time()

for random_i in range(1, Nrealisations):
    # Random borders generation
    borders_random = []
    p = 0
    for i in range(0, len(b)):
        sep = int(
            np.random.normal(
                loc=mean(list_size), scale=std(list_size), size=None
            )
        )
        p = p + sep
        if p < n1:
            borders_random.append(p)
    np.savetxt(
        "Borders_realisation" + "_" + str(random_i) + ".txt",
        borders_random,
        fmt="%d",
    )

    # Adding of CIDs/TADs borders:--------------------------------------------
    MAT_CIDs = np.ones((shape(m)))
    area = int(shape(RATIO_CIDs)[0] / 2)
    MAT_CIDs = np.concatenate((MAT_CIDs, MAT_CIDs, MAT_CIDs), axis=1)
    MAT_CIDs = np.concatenate((MAT_CIDs, MAT_CIDs, MAT_CIDs), axis=0)
    nb = 0
    for bi in borders_random:
        bi = int(bi)
        nb += 1
        MAT_CIDs_i = np.ones((shape(MAT_CIDs)))
        MAT_CIDs_i[
            np.ix_(
                range(n1 + bi - area, n1 + bi + area + 1),
                range(n1 + bi - area, n1 + bi + area + 1),
            )
        ] = RATIO_CIDs
        MAT_CIDs = MAT_CIDs * MAT_CIDs_i
    MAT_CIDs = MAT_CIDs[np.ix_(range(n1, 2 * n1), range(n1, 2 * n1))]
    MAT_CIDs = (MAT_CIDs + np.transpose(MAT_CIDs)) / 2  #  resymetrisation

    # Adding of LOOPS patterns:  ----------------------------------------------
    combi_pos = list(itertools.combinations(borders_random, 2))
    loops_random = random.sample(combi_pos, 10)

    loops_random = []
    for b1 in borders_random:
        for b2 in borders_random:
            if b2 - b1 > 2 and b2 - b1 < 100:  #  loops detectable until 200 kb
                if (100 - (b2 - b1)) / 100.0 > rand() * 5.0:
                    loops_random.append((b1, b2))
                    print(b2 - b1)
    loops_random = set(loops_random)
    np.savetxt(
        "Loops_realisation" + "_" + str(random_i) + ".txt",
        list(loops_random),
        fmt="%d",
    )
    MAT_LOOPS = np.ones((shape(m)))
    area = int(shape(RATIO_LOOPS)[0] / 2)
    MAT_LOOPS = np.concatenate((MAT_LOOPS, MAT_LOOPS, MAT_LOOPS), axis=1)
    MAT_LOOPS = np.concatenate((MAT_LOOPS, MAT_LOOPS, MAT_LOOPS), axis=0)
    for l in loops_random:
        MAT_LOOPS_i = np.ones((shape(MAT_LOOPS)))
        MAT_LOOPS_i[
            np.ix_(
                range(n1 + int(l[0]) - area, n1 + int(l[0]) + area + 1),
                range(n1 + int(l[1]) - area, n1 + int(l[1]) + area + 1),
            )
        ] = RATIO_LOOPS
        MAT_LOOPS = MAT_LOOPS * MAT_LOOPS_i
        MAT_LOOPS_i[
            np.ix_(
                range(n1 + int(l[1]) - area, n1 + int(l[1]) + area + 1),
                range(n1 + int(l[0]) - area, n1 + int(l[0]) + area + 1),
            )
        ] = np.transpose(RATIO_LOOPS)
        MAT_LOOPS = MAT_LOOPS * MAT_LOOPS_i

    MAT_LOOPS = MAT_LOOPS[
        np.ix_(range(n1, 2 * n1), range(n1, 2 * n1))
    ]  #  we refocus
    MAT_LOOPS = (MAT_LOOPS + np.transpose(MAT_LOOPS)) / 2  #  resymetrisation
    #    imshow(MAT_LOOPS, vmin=0.0, vmax = 2.0, interpolation ="none", cmap = "seismic")
    #    colorbar()

    # Complete Propensities Map with all features :
    #    MAT_PROPEN = MAT_GENO * MAT_CIDs
    #    MAT_PROPEN = MAT_GENO * MAT_LOOPS
    MAT_PROPEN = MAT_GENO * MAT_CIDs * MAT_LOOPS
    MAT_PROPEN = scn.scn_func(MAT_PROPEN, 0)

    # Conversion of the propensities matrice into a vector :
    VECT_PROPEN = np.reshape(MAT_PROPEN, (s * s, 1))
    VECT_PROPEN = VECT_PROPEN * 1.0 / sum(VECT_PROPEN)
    VECT_PROPEN = np.squeeze(VECT_PROPEN)
    VECT_INDICES = range(0, s * s)
    custm = stats.rv_discrete(name="custm", values=(VECT_INDICES, VECT_PROPEN))

    Nreads2 = int(Nreads / 2)
    VECT_REALISATION = custm.rvs(size=Nreads2)  #  time consuming !!

    # Come back to the matrice:
    MAT_SIMUL = np.zeros((shape(m)))
    for v in VECT_REALISATION:
        (i, j) = VECT_COMPT[v]
        MAT_SIMUL[i, j] += 1
        MAT_SIMUL[j, i] += 1

    np.savetxt(
        "MAT_RAW_realisation" + "_" + str(random_i) + ".txt",
        MAT_SIMUL,
        fmt="%d",
    )

    matscn2 = scn.scn_func(MAT_SIMUL, 0)
    matscn2 = (matscn2 + np.transpose(matscn2)) / 2.0
    np.savetxt(
        "MAT_NORMALISED_realisation" + "_" + str(random_i) + ".txt",
        matscn2,
        fmt="%e",
    )

    imshow(
        matscn2 ** 0.15,
        interpolation="none",
        cmap="afmhot_r",
        vmin=0.0,
        vmax=0.8,
    )
    plt.savefig(
        "MAT_realisation" + str(random_i) + ".png", dpi=500, format="png"
    )


t2 = time.time()
print(str(int(t2 - t1)) + " sec")
