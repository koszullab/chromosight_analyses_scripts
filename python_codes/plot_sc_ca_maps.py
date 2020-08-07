# Show one matrix which was used to build a chromosight kernel and circle pixels used in the kernel
# Show a second matrix on the right on which detection was performed, circle the detected pattern
# cmdoret, 20200129

import sys
import matplotlib.pyplot as plt
import cooler
import numpy as np
import pandas as pd

# load files
cool_1 = cooler.Cooler(sys.argv[1])
cool_2 = cooler.Cooler(sys.argv[2])
in_coords = pd.read_csv(sys.argv[3], sep='\t', header=None)
out_coords = pd.read_csv(sys.argv[4], sep='\t')

# Generate display
fig, ax = plt.subplots(1, 2)
ax[0].imshow(np.log(cool_1.matrix()[:]), cmap='afmhot_r')
ax[0].scatter(in_coords[1], in_coords[0])
ax[0].set_title('S. cerevisiae')
ax[1].imshow(np.log(cool_2.matrix()[:]), cmap='afmhot_r')
ax[1].scatter(out_coords.bin1, out_coords.bin2, facecolors='none', edgecolors='blue', s=20)
ax[1].set_title('C. albicans')
plt.show()

fig, ax = plt.subplots(1, in_coords.shape[0])
for i, (x, y) in in_coords.iterrows():
    ax[i].imshow(np.log(cool_1.matrix()[x-20:x+20, y-20:y+20]), cmap='afmhot_r')
plt.show()
