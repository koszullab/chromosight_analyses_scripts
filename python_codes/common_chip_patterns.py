
"""
This script takes a multiple files with 2d coordinates (2 columns). The first file is considered the
"reference" file, and all coordinate intersections between the reference and each subsequent files
are computed. A JITTER value can be specified to allow errors when intersecting coordinates. All
coordinates are first smeared by +/- JITTER pixels before computing the intersect of the resulting spots.

20200814, cmdoret
""" 
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import common_patterns_upsetplot as comm

JITTER = 1

nargs = len(sys.argv) - 1
set_names, all_patterns = comm.load_coord_files(sys.argv[1:])
ref = set_names[0]
jittered = [None] * nargs
# Find the pattern with the highest coordinate in x or y among all samples.
# This (with added jitter) will be used as the matrix shape.
mat_shape = max([max(m.iloc[:, 0].max(), m.iloc[:, 1].max()) for m in all_patterns]) + 1 + JITTER
for i, pat in enumerate(all_patterns):
    jittered[i] = comm.jitter_coords(pat, jitter=JITTER, shape=mat_shape)

inter_df = comm.multi_intersect_spots(set_names, jittered)

# Only keep sets which include the reference sample
ref_df = inter_df.loc[inter_df[ref], :].drop(columns=[ref])
degree = inter_df.iloc[:, :-1].sum(axis=1)
# Only keep 2-way intersections with ref 
two_way_df = ref_df.loc[degree == 2, :]
# Keep self counts (no intersections)
one_way_df = inter_df.loc[(degree == 1) & (~ inter_df[ref]), :].drop(columns=[ref])
two_way_df['proportion'] = 0

# Compute the proportion of patterns overlapping ref, for each software
for i in set_names[1:]:
    n_inter_ref = two_way_df.loc[two_way_df[i], 'n_patterns'].values[0]
    n_total = one_way_df.loc[one_way_df[i], 'n_patterns'].values[0] 
    two_way_df.loc[two_way_df[i], 'proportion'] =  n_inter_ref / n_total

# Reformat dataframe to long format for plotting
plot_df = two_way_df.melt(id_vars=['n_patterns', 'proportion'], var_name='software')
plot_df = plot_df.loc[plot_df['value']].drop(columns=['value'])

plt.bar(x=plot_df.software, height=plot_df.proportion*100)
plt.title("Overlap of loop calls with CTCF peaks")
plt.ylabel("Percent loop overlapping CTCF peaks")
plt.show()
