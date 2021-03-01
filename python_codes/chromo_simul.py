# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:27:26 2018
7: reprise feb 2020
Notably we try other pattern of borders (the generic ones) to make propensity maps.
"""
import random
import itertools
import os
from os.path import join
import sys
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import cooler
import hicstuff.hicstuff as hcs
import chromosight.utils.preprocessing as cup
import chromosight.kernels as ck

# Parse CLI arguments
# Path to the cool file with experimental data
clr = cooler.Cooler(sys.argv[1])
# Read selected chromosome
chrom = sys.argv[2]
# Path to positions of borders in the experimental matrix (in bins,
# relative to chromosome start)
borders_pos = np.loadtxt(sys.argv[3])
# Create output directory if it does not exist already
out_dir = Path(sys.argv[4])
out_dir.mkdir(parents=True, exist_ok=True)

# Load experimental (balanced !) Hi-C data:
# To balance it, you can use cooler's CLI command on the file: cooler balance
mat = clr.matrix(balance=True, sparse=False).fetch(chrom)
# Load genomic coordinates of bins (make indices relative to chromosome)
bins = clr.bins().fetch(chrom).reset_index(drop=True)
# Retrieve missing bins that were excluded from the normalisation
bad_bins = bins.loc[np.isnan(bins.weight)].index.values
# Number of raw contacts in raw the matrix
nreads = clr.matrix(sparse=True, balance=False).fetch(chrom).sum()
# Desired number of simulated matrices to generate
Nrealisations = 10
n = mat.shape[0]  # size of the chromosome contact map

# Load patterns
# Generic pattern:
ratio_borders = ck.borders["kernels"][0]
ratio_borders = ratio_borders ** 0.2  # to attenuate a bit
ratio_loops = ck.loops["kernels"][0]
ratio_loops = ratio_loops  # to attenuate a bit
# Loops and borders windows radii
l_rad = int(ratio_loops.shape[0] // 2)
b_rad = int(ratio_borders.shape[0] // 2)

# Remove borders from experimental data overlapping empty bins:
borders_pos = list(borders_pos)
for b in borders_pos.copy():
    # If any bad (empty) bin is closer to a border than the radius of the
    # pattern template, drop it
    if np.any(np.abs(bad_bins - b) < max(ratio_borders.shape)):
        borders_pos.remove(b)


def vec1d_to_2d(pos, n_cols):
    """
    Given a 1-dimensional position, return a tuple (row,col) indicating the
    corresponding coordinate in a matrix of width n_cols.
    """
    return pos // n_cols, pos % n_cols


def vec2d_to_1d(row, col, n_cols):
    """
    Given row and column indices for a matrix of width n_cols, return the
    corresponding 1-dimensional index.
    """
    return row * n_cols + col


#  Genomic distance law from experimental data:
d_smooth = cup.distance_law(
    clr.matrix(sparse=True, balance=True).fetch(chrom).tocsr(),
    detectable_bins=~np.isnan(bins.weight),
    smooth=False,
)
d_smooth[np.isnan(d_smooth)] = 0
prob_d_smooth = d_smooth / np.sum(d_smooth)
# TADs size distribution:
borders_sizes = []
for i in range(1, len(borders_pos)):
    bi = borders_pos[i] - borders_pos[i - 1]
    borders_sizes.append(bi)

# Genomic distance law:
mat_geno = np.zeros(mat.shape)
for i in range(0, n):
    for j in range(0, n):
        mat_geno[i, j] = prob_d_smooth[abs(j - i)]

mat_geno[mat_geno < 0] = 0


# ------------------------------------------------------------------------------
#  RANDOM MAPs GENERATION
# ------------------------------------------------------------------------------
# NOTE: Send to output folder
t1 = time.time()

for random_i in range(1, Nrealisations):
    # Random borders generation
    borders_random = []
    p = 0
    for i in range(0, len(borders_pos)):
        sep = int(
            np.random.normal(
                loc=np.mean(borders_sizes),
                scale=np.std(borders_sizes),
                size=None,
            )
        )
        p += sep
        if 0 < p < n:
            borders_random.append(p)
    np.savetxt(
        join(out_dir, f"Borders_realisation_{random_i}.txt"),
        borders_random,
        fmt="%d",
    )

    # Adding of CIDs/TADs borders:--------------------------------------------
    mat_borders = np.ones(mat.shape)
    nb = 0
    for bi in borders_random:
        bi = int(bi)
        nb += 1
        mat_borders_i = np.ones(mat_borders.shape)
        try:
            mat_borders_i[
                np.ix_(
                    range(bi - b_rad, bi + b_rad + 1),
                    range(bi - b_rad, bi + b_rad + 1),
                )
            ] = ratio_borders
        # Skip if the border window goes out of the map
        except IndexError:
            continue
        mat_borders *= mat_borders_i
    # Make the matrix symmetric
    mat_borders = (mat_borders + np.transpose(mat_borders)) / 2

    # Adding of loops patterns:  ----------------------------------------------
    combi_pos = list(itertools.combinations(borders_random, 2))
    loops_random = random.sample(combi_pos, min(10, len(combi_pos)))

    # Generate loops coordinates aligned with borders
    loops_random = []
    for b1 in borders_random:
        for b2 in borders_random:
            lsize = b2 - b1
            #  loops detectable between 4 and 200 kb
            if 4000 / clr.binsize > lsize < 200000 / clr.binsize:
                # Give a probability of using each loop inversely proportional
                # to the length of the loop (-> more small loops)
                if (100 - lsize) / 100.0 > np.random.rand() * 5.0:
                    loops_random.append((b1, b2))
    loops_random = set(loops_random)
    np.savetxt(
        join(out_dir, f"Loops_realisation_{random_i}.txt"),
        list(loops_random),
        fmt="%d",
    )

    mat_loops = np.ones(mat.shape)
    for loop in loops_random:
        loop = np.array(loop).astype(int)
        mat_loops_i = np.ones(mat_loops.shape)
        # Note we fill the upper triangle (l[0] > l[1])
        try:
            mat_loops_i[
                np.ix_(
                    range(loop[1] - l_rad, loop[1] + l_rad + 1),
                    range(loop[0] - l_rad, loop[0] + l_rad + 1),
                )
            ] = ratio_loops
        # Skip if either loop anchor goes out of the map
        except IndexError:
            continue
        mat_loops *= mat_loops_i
    # Set lower triangle to 0
    mat_loops = np.triu(mat_loops)
    # Make matrix symmetric by adding upper triangle to lower triangle
    mat_loops = mat_loops + np.transpose(mat_loops)
    # Remove half of diagonal since it has been doubled during transpose
    mat_loops[np.diag_indices_from(mat_loops)] /= 2

    # The final propensity map will be the (elementwise) product of the genomic
    # distance law map, the border map and the loop map:
    mat_propen = mat_geno * mat_borders * mat_loops
    mat_propen = hcs.normalize_dense(mat_propen, iterations=40)
    # Conversion of the propensities matrice into a 1D vector :
    vect_propen = np.reshape(mat_propen, (n * n, 1))
    # Convert to probabilities
    vect_propen = vect_propen * 1.0 / sum(vect_propen)
    vect_propen = np.squeeze(vect_propen)
    vect_indices = range(0, n * n)
    # Generate a distribution from the vector  of probabilities and sample contacts

    # Use the probability distribution to sample contacts, time consuming !!
    mat_simul = np.zeros(mat.shape)
    nreads_left = nreads
    while nreads_left:
        # Fill by batches of 100k contacts to spare memory
        sample_batch = np.random.choice(
            vect_indices, size=min(100000, nreads_left), p=vect_propen
        )
        batch_rows, batch_cols = vec1d_to_2d(sample_batch, n)
        mat_simul[batch_rows, batch_cols] += 1
        print(f"{nreads - nreads_left} / {nreads} contacts subsampled")
        nreads_left -= len(sample_batch)

    # Symmetrise matrix
    mat_simul += np.transpose(mat_simul)
    mat_simul /= 2

    np.savetxt(
        join(out_dir, f"MAT_RAW_realisation_{random_i}.txt"),
        mat_simul,
        fmt="%d",
    )

    matscn2 = hcs.normalize_dense(mat_simul, iterations=60)
    np.savetxt(
        join(out_dir, f"MAT_NORMALISED_realisation{random_i}.txt"),
        matscn2,
        fmt="%e",
    )
    plt.imshow(
        matscn2,
        interpolation="none",
        cmap="afmhot_r",
    )
    plt.savefig(
        join(out_dir, f"MAT_realisation{random_i}.png"), dpi=500, format="png",
    )


t2 = time.time()
print(str(int(t2 - t1)) + " sec")
