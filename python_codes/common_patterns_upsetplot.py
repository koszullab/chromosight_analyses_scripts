"""
This script taks multiple files with loop coordinates (2 columns) and computes
The intersect between each combination of set.
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

jitter = 1
nargs = len(sys.argv) - 1
all_patterns = [None] * nargs
set_names = [None] * nargs
jittered = [None] * nargs

for i, in_file in enumerate(sys.argv[1:]):    
    # Load coordinates into dataframes
    all_patterns[i] = pd.read_csv(in_file, header=None, names=['bin1', 'bin2'], sep='\t')
    # Remember set name based on filename
    set_names[i] = os.path.basename(in_file)

# Find the pattern with the highest coordinate in x or y among all samples.
# This (with added jitter) will be used as the matrix shape.
mat_shape = max([max(m.iloc[:, 0].max(), m.iloc[:, 1].max()) for m in all_patterns]) + 1 + jitter

for i, pat in enumerate(all_patterns):
    n_coords = pat.shape[0]
    # Allocate an empty vector to store all jittered coordinates
    x_jitter = np.zeros(n_coords * (2*jitter+1)**2, dtype=int)
    y_jitter = np.zeros(n_coords * (2*jitter+1)**2, dtype=int)
    
    start = 0
    # Compute coordinates with each jitter combination
    for x in range(-jitter, jitter+1):
        for y in range(-jitter, jitter+1):
            x_jitter[start:start+n_coords] = pat.iloc[:, 0] + x
            y_jitter[start:start+n_coords] = pat.iloc[:, 1] + y
            start += n_coords
    # Generate a sparse matrix with jittered 2D coords
    jittered[i] = sp.coo_matrix(
        (np.ones(len(x_jitter)), (x_jitter, y_jitter)),
        shape=(mat_shape, mat_shape),
        dtype=int
    )
    # Make sure values are either 0 (absent) or 1 (present)
    jittered[i].data[jittered[i].data > 1] = 1
    
# Generate all set combinations
def combos(size): return [i for i in it.combinations(range(nargs), size)]
all_combos = []

for size in range(1, nargs + 1):
    all_combos += combos(size)

n_combos = len(all_combos)

set_df = pd.DataFrame({set_names[i]: np.zeros(n_combos, dtype=bool) for i in range(nargs)})
set_df['n_patterns'] = 0

# Loop over each combo
for c, combo in enumerate(all_combos):
    # intersection (AND) of all matrix elements in combo
    for e, elem in enumerate(combo):
        # Use pixel-wise minimum to compute intersection
        if not e:
            inter = jittered[elem]
        else:
            inter = inter.minimum(jittered[elem])
    # Retrieve original pattern coordinates matching intersection
    num_spots, _ = cud.label_foci(inter)
    # Store number of common patterns (intersection spots) for current combo in dataframe
    set_df.loc[c, [set_names[i] for i in combo]] = True
    set_df.loc[c, 'n_patterns'] = num_spots

print(set_df)
# Reformat dataframe into a series with nested indices, for upsetplot com   patibility
sets = set_df.set_index([set_names[i] for i in range(nargs)]).n_patterns
upsetplot.plot(sets, show_counts=True)
plt.show()
