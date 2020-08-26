# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 14:27:26 2018
7: reprise feb 2020
Notably we try other pattern of borders (the generic ones) to make propensity maps.
"""
import scipy.ndimage
import scipy.io as sio
from pylab import *
import random
import itertools

import scipy.stats as stats
from scipy.stats import rv_discrete

import glob
import os
import time

import statsmodels.api as sm
lowess = sm.nonparametric.lowess
import matplotlib

os.chdir("/home/axel/Bureau/z_python_scripts_copy/")
import scn
import ice_mirny3
import distance_law_human
import hicstuff as hcs

# Experimental data:
m=loadtxt("/home/axel/Bureau/YEAST/data_yeast_pds5/ALL_EXPERIMENTS_PDF/AT147_Pds5-AID-noTir1-G1-cdc20-TS_files/txt/MAT_RAW_chr5_AT147_Pds5-AID-noTir1-G1-cdc20-TS_2kb.txt")
Nreads = sum(m)
Nrealisations = 10

#imshow(m**0.15,interpolation="none",cmap= 'afmhot_r')
#colorbar()
n1 = shape(m)[0]   # size of the chromosome contact map 

# Normalisation:
matscn1=ice_mirny3.ice_func(m,0)
matscn1=scn.scn_func(matscn1,0)


# Borders from experimental data:
chrm=5
b = loadtxt("/home/axel/Bureau/YEAST/data_yeast_pds5/borders_detection/borders_POISSIAN_Dplus_chr"+str(chrm)+"_.txt", dtype="int")
b = list(b)
ele_to_remove = [221,224,246,249]
for e in ele_to_remove :
    if e in b:
        b.remove(e)

#plot(b, b, 'o', color="cyan")

# To convert the 1D vetor into matrice object:
MAT_INDICES = np.zeros(( shape(m)  ))
compt = 0
VECT_COMPT = []
for i, j in itertools.product(range(n1), range(n1)):
    MAT_INDICES[i,j] = compt
    VECT_COMPT.append( (i,j) )
    compt += 1
#imshow(MAT_INDICES, interpolation ="none")


#  Genomic distance law from experiemental data: 
s = shape(matscn1)[0]
print(s)
prob_d = distance_law_human.dist_law(matscn1)
prob_d =  prob_d / sum(prob_d)
#plot(prob_d)

# LOWESS Fit:
z = lowess(prob_d[:,0], range(0, len(prob_d) ) , frac=1./15)
plot(z[:,0],z[:,1],lw=3)
prob_d_lowess = z[:,1]

# TADs size distribution:
list_size=[]
for i in range(1,len(b) ) :
    bi = b[i] - b[i-1]
    list_size.append(bi)
mean(list_size) 
std(list_size) 
min(list_size)
max(list_size) 

s = shape(matscn1)[0]

# Genomic distance law:
MAT_GENO = np.zeros(( shape(m)  ))
for i in range(0,s) :
    for j in range(0,s) :
        MAT_GENO[i,j] = prob_d_lowess[ abs(j-i) ]
 
if(len(MAT_GENO[MAT_GENO<0]) >0) :
    print("attention presence of neg elements")
    print(len(MAT_GENO[MAT_GENO<0]) )  
MAT_GENO[MAT_GENO<0]=0       

# Generic pattern:
RATIO_CIDs = loadtxt("/home/axel/Bureau/tools/chromosight/chromosight/kernels/artificial_template_borders_type1.txt")
RATIO_CIDs = RATIO_CIDs**0.2  # to attenuate a bit
imshow(RATIO_CIDs, vmin= 0.0, vmax = 2.0, interpolation ="none", cmap = "seismic")
colorbar()
shape(RATIO_CIDs)

RATIO_LOOPS = loadtxt("/home/axel/Bureau/chromovision-master/data/artificial_template_loops_type1.txt") # too strong 
RATIO_LOOPS = RATIO_LOOPS**0.7  # to attenuate a bit
imshow(RATIO_LOOPS, vmin=0.0, vmax = 2.0, interpolation ="none", cmap = "seismic")
colorbar()

close("all")

#------------------------------------------------------------------------------
#  RANDOM MAPs GENERATION 
#------------------------------------------------------------------------------
os.chdir("/home/axel/Bureau/TRAINING_SET_15_2")
t1 = time.time()

for random_i in range(1, Nrealisations) :
    # Random borders generation 
    borders_random = []
    p = 0
    for i in range(0,len(b) ) :
        sep = int( np.random.normal(loc=mean(list_size), scale=std(list_size), size=None ) )
        p = p + sep
        if p < n1 :
            borders_random.append(p)
    len(borders_random)
    np.savetxt("Borders_realisation"+"_"+str( random_i )+".txt", borders_random, fmt='%d')
    
    # Adding of CIDs/TADs borders:--------------------------------------------
    MAT_CIDs = np.ones(( shape(m)  ))
    shape(MAT_CIDs)
    area = int(shape(RATIO_CIDs)[0] / 2 )
    MAT_CIDs = np.concatenate((MAT_CIDs, MAT_CIDs, MAT_CIDs), axis=1)
    MAT_CIDs = np.concatenate((MAT_CIDs, MAT_CIDs, MAT_CIDs), axis=0)
    shape(MAT_CIDs)
    nb = 0
    for bi in borders_random :
        print(bi)
        bi=int(bi)
        nb += 1
        MAT_CIDs_i = np.ones(( shape(MAT_CIDs)  ))
        MAT_CIDs_i[ np.ix_(range(n1+bi-area,n1+bi+area+1)  , range(n1+bi-area,n1+bi+area+1) ) ] = RATIO_CIDs
        MAT_CIDs = MAT_CIDs * MAT_CIDs_i
    #MAT_CIDs = MAT_CIDs / nb
    MAT_CIDs = MAT_CIDs[ np.ix_(range(n1,2*n1)  , range(n1,2*n1)   ) ]
    MAT_CIDs = (MAT_CIDs  + np.transpose(MAT_CIDs ) ) / 2  #  resymetrisation
    shape(MAT_CIDs)
#    imshow(MAT_CIDs, vmin=0.0, vmax = 2.0, interpolation ="none", cmap = "seismic")
#    colorbar()
    imshow(MAT_CIDs**0.2,  interpolation ="none", cmap = "afmhot_r")
    imshow(MAT_CIDs,  interpolation ="none", cmap = "seismic", vmin=0.0, vmax =2.0)
    
    # Adding of LOOPS patterns:  ----------------------------------------------
    combi_pos = list(itertools.combinations(borders_random,2) )
    len(combi_pos)    
    loops_random = random.sample(combi_pos , 10)
    
    loops_random = []
    for b1 in borders_random :
        for b2 in borders_random :
            if b2 - b1 > 2 and b2 - b1 < 100 :   #  loops detectable until 200 kb
                if ( 100 - (b2-b1) ) / 100.0  >  rand() * 5. :
                    loops_random.append( (b1, b2) )
                    print(b2-b1)
    len(loops_random)
    loops_random = set(loops_random)
    len(loops_random)
    np.savetxt("Loops_realisation"+"_"+str( random_i )+".txt", list(loops_random) , fmt='%d')
    MAT_LOOPS = np.ones(( shape(m)  ))
    area = int(shape(RATIO_LOOPS)[0] / 2)
    MAT_LOOPS = np.concatenate((MAT_LOOPS, MAT_LOOPS, MAT_LOOPS), axis=1)
    MAT_LOOPS = np.concatenate((MAT_LOOPS, MAT_LOOPS, MAT_LOOPS), axis=0)
    shape(MAT_LOOPS)
    for l in loops_random :
        print(l)
        MAT_LOOPS_i = np.ones(( shape(MAT_LOOPS)  ))
        MAT_LOOPS_i[ np.ix_(range(n1+int(l[0])-area, n1+int(l[0])+area+1)  , range(n1+int(l[1])-area, n1+int(l[1])+area+1) ) ] = RATIO_LOOPS
        MAT_LOOPS = MAT_LOOPS * MAT_LOOPS_i
        MAT_LOOPS_i[ np.ix_(range(n1+int(l[1])-area, n1+int(l[1])+area+1)  , range(n1+int(l[0])-area, n1+int(l[0])+area+1) ) ] = np.transpose(RATIO_LOOPS)
        MAT_LOOPS = MAT_LOOPS * MAT_LOOPS_i
        
    MAT_LOOPS = MAT_LOOPS[ np.ix_(range(n1,2*n1)  , range(n1,2*n1)   ) ]  #  we refocus 
    MAT_LOOPS = (MAT_LOOPS + np.transpose(MAT_LOOPS) ) / 2  #  resymetrisation
#    imshow(MAT_LOOPS, vmin=0.0, vmax = 2.0, interpolation ="none", cmap = "seismic")
#    colorbar()
    
    # Complete Propensities Map with all features :
#    MAT_PROPEN = MAT_GENO * MAT_CIDs
#    MAT_PROPEN = MAT_GENO * MAT_LOOPS
    MAT_PROPEN = MAT_GENO * MAT_CIDs * MAT_LOOPS
    MAT_PROPEN = scn.scn_func(MAT_PROPEN, 0)
    #plot(MAT_PROPEN.sum(axis=0) )
    imshow(MAT_PROPEN**.25, interpolation = "none" ,cmap= 'jet', vmin=0.0, vmax= 0.8 )
    colorbar()
    imshow(MAT_PROPEN**.2, interpolation = "none" ,cmap= 'afmhot_r', vmin=0.0, vmax= 0.8 )
    colorbar()
    for l in loops_random :
        plt.scatter(l[0], l[1], s=80, facecolors='none', edgecolors='yellow')
    
    # Conversion of the propensities matrice into a vector :
    VECT_PROPEN = np.reshape(MAT_PROPEN, (s*s, 1))
    VECT_PROPEN = VECT_PROPEN * 1. / sum(VECT_PROPEN)
    VECT_PROPEN=np.squeeze(VECT_PROPEN)
    sum(VECT_PROPEN)
    VECT_INDICES = range(0, s*s )
    custm  = stats.rv_discrete(name='custm', values=(VECT_INDICES, VECT_PROPEN))
    
    Nreads2 = int( Nreads/ 2) 
    VECT_REALISATION = custm.rvs(size=Nreads2)    #  time consuming !!
    
    # Come back to the matrice:
    MAT_SIMUL = np.zeros(( shape(m)  ))
    for v in VECT_REALISATION : 
        (i,j) = VECT_COMPT[v]
        MAT_SIMUL[i,j] += 1
        MAT_SIMUL[j,i] += 1
        
#    # Plots: 
#    imshow(MAT_SIMUL**0.15, interpolation = "none" ,cmap= 'afmhot_r')
#    sum(MAT_SIMUL)    
    np.savetxt("MAT_RAW_realisation"+"_"+str( random_i )+".txt", MAT_SIMUL, fmt='%d')
      
    matscn2 = scn.scn_func(MAT_SIMUL, 0)
    matscn2= (matscn2 + np.transpose(matscn2) ) / 2.0 
    np.savetxt("MAT_NORMALISED_realisation"+"_"+str( random_i )+".txt", matscn2, fmt='%e')
   
#    plot(borders_random, borders_random, 'o', color="cyan")
#    imshow(matscn2**0.15,interpolation="none",cmap= 'afmhot_r', vmin=0.0, vmax = 0.8)
#    colorbar()
#    title("Simulated Contact Map no"+str( random_i ) + " Number of reads=" + str(Nreads2) )
#    #plot(borders_random, borders_random,'o',color="cyan")
#    for l in loops_random :
#        plt.scatter(l[0], l[1], s=80, facecolors='none', edgecolors='yellow')
        
    imshow(matscn2**0.15,interpolation="none",cmap= 'afmhot_r', vmin=0.0, vmax = 0.8)    
    plt.savefig("MAT_realisation"+str( random_i )+".png", dpi=500, format='png')
    plt.close()


t2 = time.time()
print( str( int(t2-t1) ) + " sec")
#------------------------------------------------------------------------------ 
#  Cmp genomic distance laws: 
#prob_d = distance_law_human.dist_law(matscn1)
#prob_d =  prob_d / sum(prob_d)
#plot(prob_d, label="EXP")
#prob_d2 = distance_law_human.dist_law(matscn2)
#prob_d2 =  prob_d2 / sum(prob_d2)
#plot(prob_d2, label="SIMUL")
#legend()


m=np.loadtxt("/home/axel/Bureau/TRAINING_SET_9/MAT_NORMALISED_realisation_9.txt")

plt.figure(1)
imshow(m**0.15, interpolation = 'None', vmin=0.0, vmax=0.8, cmap="afmhot_r")
plt.colorbar()

plt.figure(2)
imshow(matscn1**0.15, interpolation = 'None', vmin=0.0, vmax=0.8, cmap="afmhot_r")
plt.colorbar()








