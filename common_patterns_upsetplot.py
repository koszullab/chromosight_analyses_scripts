"""
This script takes a multiple files with 2d coordinates (2 columns) and computes all coordinate
intersections between each combination of files. A JITTER value can be specified to allow errors
when intersecting coordinates. All coordinates are first smeared by +/- JITTER pixels before
computing the intersect of the resulting spots.
20200810, cmdoret
""" 
import sys
import os
import itertools as it
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import chromosight.utils.detection as cud
import upsetplot

JITTER = 1


def load_coord_files(paths):
    """
    Given a list of paths to 2-column tab-separated
    files containing 2D coordinates, return a list
    of their filenames, and a list of dataframes,
    each containing the 2D coords of a file, with
    column names ['bin1', 'bin2'].
    """
    n = len(paths)
    all_patterns = [None] * n
    set_names = [None] * n
    for i, in_file in enumerate(paths):    
        # Load coordinates into dataframes
        all_patterns[i] = pd.read_csv(in_file, header=None, names=['bin1', 'bin2'], sep='\t')
        # Remember set name based on filename
        set_names[i] = os.path.basename(in_file)
    return set_names, all_patterns


def jitter_coords(df, jitter=JITTER, shape=None):
    """
    Given a dataframes of 2 columns names ['bin1', 'bin2'],
    generate a sparse matrix containing explicit values (ones)
    at the coords stored in the input dataframe, and all their
    neighbours within a range given by jitter.
    A shape for the output matrix can be predefined, otherwise
    it is determined by the maximum coordinate and jitter value.

    Note: Having a jitter of 1 means each coordinate will be smeared
    on +/- 1 pixel around, resulting in a 3x3 pixel square.
    """
    n_coords = df.shape[0]
    # Allocate an empty vector to store all jittered coordinates
    x_jitter = np.zeros(n_coords * (2*jitter+1)**2, dtype=int)
    y_jitter = np.zeros(n_coords * (2*jitter+1)**2, dtype=int)
    
    start = 0
    # Compute coordinates with each jitter combination
    for x in range(-jitter, jitter+1):
        for y in range(-jitter, jitter+1):
            x_jitter[start:start+n_coords] = df.iloc[:, 0] + x
            y_jitter[start:start+n_coords] = df.iloc[:, 1] + y
            start += n_coords
    # Generate a sparse matrix with jittered 2D coords
    jit_mat = sp.coo_matrix(
        (np.ones(len(x_jitter)), (x_jitter, y_jitter)),
        shape=(shape, shape),
        dtype=int
    )
    # Make sure values are either 0 (absent) or 1 (present)
    jit_mat.data[jit_mat.data > 1] = 1
    return jit_mat


def multi_intersect_spots(names, mats):
    """
    Given a list of sample names and another list of their
    corresponding sparse jittered coordinate matrices, return
    a dataframe with the number of common jittered coordinates
    between each combination of samples, with degree 1 - n_samples. 
    
    Note: The input matrices should contain square "spots". Each of
    these spots correspond to 1 original jittered cooridnate. We compute
    the intersect of the number of individual spots between samples to
    allow a small error.
    """
    if len(names) != len(mats):
        raise ValueError("Number of sample names and matrices differ.")

    # Generate all set combinations
    n_samples = len(names)
    def combos(size): return [i for i in it.combinations(range(n_samples), size)]
    all_combos = []
    for size in range(1, n_samples + 1):
        all_combos += combos(size)

    n_combos = len(all_combos)

    set_df = pd.DataFrame({names[i]: np.zeros(n_combos, dtype=bool) for i in range(n_samples)})
    set_df['n_patterns'] = 0

    # Loop over each combo
    for c, combo in enumerate(all_combos):
        # intersection (AND) of all matrix elements in combo
        for e, elem in enumerate(combo):
            # Use pixel-wise minimum to compute intersection
            if not e:
                inter = mats[elem]
            else:
                inter = inter.minimum(mats[elem])
        # Retrieve original pattern coordinates matching intersection
        num_spots, _ = cud.label_foci(inter)
        # Store number of common patterns (intersection spots) for current combo in dataframe
        set_df.loc[c, [names[i] for i in combo]] = True
        set_df.loc[c, 'n_patterns'] = num_spots
    return set_df

if __name__ == '__main__':
    nargs = len(sys.argv) - 1
    set_names, all_patterns = load_coord_files(sys.argv[1:])
    jittered = [None] * nargs
    # Find the pattern with the highest coordinate in x or y among all samples.
    # This (with added jitter) will be used as the matrix shape.
    mat_shape = max([max(m.iloc[:, 0].max(), m.iloc[:, 1].max()) for m in all_patterns]) + 1 + JITTER
    for i, pat in enumerate(all_patterns):
        jittered[i] = jitter_coords(pat, jitter=JITTER, shape=mat_shape)

    set_df = multi_intersect_spots(set_names, jittered)

    print(set_df)
    # Reformat dataframe into a series with nested indices, for upsetplot com   patibility
    sets = set_df.set_index([set_names[i] for i in range(nargs)]).n_patterns
    upsetplot.plot(sets, show_counts=True)
    plt.show()
