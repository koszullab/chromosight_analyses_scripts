# -*- coding: utf-8 -*-
"""
@author: axel KournaK 
To vizualise detected loops! 
"""
import time
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import random as rando
import sys
import os
sys.path.insert(0, os.path.abspath("/home/axel/Bureau/z_python_scripts_copy"))
import scn 
import ice_mirny3
import scipy
import scipy.ndimage
import scipy.io as sio
import distance_law_human
import hicstuff as hcs
import numpy as np
import json

import sys
import chromosight.utils.plotting as cup 
import pandas as pd
import cooler

cool_file = sys.argv[1]
loop_file = sys.argv[2]
name_bank = sys.argv[3]

cool_file="/home/axel/Bureau/Arima/contacts2_1_SRR6675327.cool"
loop_file='/home/axel/Bureau/Arima/out_contacts2_1_SRR6675327.cool.tsv'
loop_file='/home/axel/Bureau/Arima/out_2_contacts2_1_SRR6675327.cool.tsv'
name_bank= "Arima kernel1"

cool_file="/media/axel/RSG4/in_situ_ChIA-PET_CTCF/4DNFIMH3J7RW.mcool.5000"
loop_file='/media/axel/RSG4/in_situ_ChIA-PET_CTCF/out_4DNFIMH3J7RW.mcool.5000.tsv'
loop_file='/media/axel/RSG4/in_situ_ChIA-PET_CTCF/out_loops1_4DNFIMH3J7RW.mcool.5000.tsv'
loop_file='/media/axel/RSG4/in_situ_ChIA-PET_CTCF/out_pear05_iterations_4DNFIMH3J7RW.mcool.5000.tsv'
loop_file='/media/axel/RSG4/chia_drop_GM/out_2_4DNFITPJZ3AW.mcool.5000.tsv'
name_bank= "in_situ_ChIA-PET_CTCF_5000_kernel1_chia-pet"

cool_file="/media/axel/RSG4/chia_drop_GM/4DNFI81RQ431.mcool.5000"
cool_file="/media/axel/RSG4/chia_drop_GM/4DNFITPJZ3AW.mcool.5000"
loop_file='/media/axel/RSG4/chia_drop_GM/out_pear03_4DNFI81RQ431.mcool.5000.tsv'
loop_file='/media/axel/RSG4/chia_drop_GM/out_2_4DNFITPJZ3AW.mcool.5000.tsv'
name_bank= "ChIA-Drop small loops_5000 kernel1"

cool_file="/media/axel/RSG4/sprite_GM12878/4DNFIUOOYQC3_gm12878_sprite.mcool.5000"
loop_file='/media/axel/RSG4/sprite_GM12878/out_4DNFIUOOYQC3_gm12878_sprite.mcool.10000.tsv'
loop_file='/media/axel/RSG4/sprite_GM12878/out_loops1_4DNFIUOOYQC3_gm12878_sprite.mcool.5000.tsv'
name_bank= "SPRITE kernel1"

cool_file="/media/axel/RSG4/HiChip_cohesin_GM/GSE80820_HiChIP_GM_cohesin_10000.cool"
loop_file='/media/axel/RSG4/HiChip_cohesin_GM/out_GSE80820_HiChIP_GM_cohesin_10000.cool.tsv'
loop_file='/media/axel/RSG4/HiChip_cohesin_GM/out_loops1_2_GSE80820_HiChIP_GM_cohesin_10000.cool.tsv'
name_bank= "HiChIP kernel1"

cool_file="/media/axel/RSG4/microC_human/4DNFI9FVHJZQ.mcool.5000"
loop_file="/media/axel/RSG4/microC_human/out_4DNFI9FVHJZQ.mcool.5000.tsv"
name_bank= "MicroC_HFFc6_expo015"

#
cool_file="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/micro-C_GM12878_redone/contacts2_micro-C_GM12878_MQ30_10000_df.cool"
loop_file="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/micro-C_GM12878_redone/out_contacts2_micro-C_GM12878_10000.tsv"
name_bank= "micro-C_GM12878_dovetail"

cool_file="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/micro-C_GM12878_redone/OmniC_800M_10000.hic.cool"
loop_file="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/micro-C_GM12878_redone/out_OmniC_800M_10000.hic.tsv"
name_bank= "OmniC_GM12878_dovetail"

cool_file="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/micro-C_GM12878_redone/CTCF.mcool.8000.cool"
loop_file="/media/axel/d0a28364-6c64-4f8e-9efc-f332d9a0f1a9/micro-C_GM12878_redone/out_CTCF.mcool.8000.cool.tsv"
name_bank= "Dovetail_CTCF.mcool.8000"

c = cooler.Cooler(cool_file)
cooler.balance_cooler(c, store=True)  # normalisation

# Input of loops detected by chromosight: 
bin_matrice = 8000 # size of bin of matrice (in bp)

group2=pd.read_table(loop_file, delimiter="\t") 
group2['start1'] = [ int(x/bin_matrice) for x in group2['start1']]
group2['start2'] = [ int(x/bin_matrice) for x in group2['start2']]
print(len(group2))

# Positions of genes of interest: 
geno_pos = pd.read_table('/home/axel/Bureau/test_chromosight/chromo_sparse_old_stuffs/interesting_human_genes.txt'
                         ,header=None, delimiter=" ")

geno_pos = pd.read_table('/home/axel/Bureau/test_chromosight/chromo_sparse_old_stuffs/interesting_human_genes.txt.lifted'
                         ,header=None, delimiter="\t")

#geno_pos = pd.read_table('/home/axel/Bureau/genes_interesting_hg38.txt2'
#                         ,header=None, delimiter=" ")

area = 1000000   #  area arount the gene of interest 
bin_size = bin_matrice # in bp 

for i in range(len(geno_pos) ) :
    print(geno_pos[0][i])
    chrm = str(geno_pos[0][i])
    pos1 = geno_pos[1][i]
    pos2 = geno_pos[2][i]
    gene = geno_pos[3][i]
    
    mat = c.matrix().fetch(region=str(chrm)+"")  # extraction of the mat from cool file
    
#    mat= c.matrix().fetch(region="chr2:190,000,000-200,000,000")
    
#    cup.plot_whole_matrix(c, loops, region="2:190,000,000-200,000,000")
    mat[np.isnan(mat)] = 0.0
    mat[np.isinf(mat)] = 0.0
    
    p1 = int( (pos1-area) / bin_size)
    p2 = int( (pos2+area) / bin_size)
    mat = mat[p1:p2,p1:p2]
    chrom_pattern_coords = group2.loc[(group2['chrom1'] == chrm )]
    chrom_pattern_coords2 =chrom_pattern_coords[ ((chrom_pattern_coords['start1'] > p1) & (chrom_pattern_coords['start1'] < p2) ) & 
                                                ((chrom_pattern_coords['start2'] > p1) & (chrom_pattern_coords['start2'] < p2) )]
    print("Number of loops in submatrice")
    print( len(chrom_pattern_coords2) )
    
    plt.imshow(np.power(mat,0.15) , interpolation="none", 
               extent=[pos1-area,pos2+area,pos2+area,pos1-area], cmap="afmhot_r")
    plt.title(gene+" - "+str(chrm)+" "+str(pos1)+" BIN="+str(bin_size) )
    
    plt.xlabel("Position along the chromosome")
    plt.colorbar()
    plt.scatter((chrom_pattern_coords2['start1']+1)*bin_size, (chrom_pattern_coords2['start2']+1)*bin_size,
                s=20, facecolors='none', edgecolors='blue')
    
    plt.savefig(gene+"_"+str(bin_size)+name_bank+".pdf")
    plt.close()



