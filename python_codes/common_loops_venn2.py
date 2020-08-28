#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 11:10:32 2019
@author: axel KournaK
To find common loops from detected coordinates with some jitter
"""
import numpy as np 
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import sys
import chromosight.utils.preprocessing as preproc

# sparse convolution
def xcorr2_sparse(signal, kernel, threshold=1e-6) :
    """
    Cross correlate a sparse 2D signal with a dense 2D kernel.
    Parameters
    ----------
    signal: scipy.sparse.csr_matrix
        A 2-dimensional numpy array Ms x Ns acting as the detrended Hi-C map.
    kernel: numpy.array of floats or tuple of numpy.arrays
        A 2-dimensional numpy array Mk x Nk acting as the pattern template. Can
        also be a factorised kernel.
    threshold : float
        Convolution score below which pixels will be set back to zero to save
        on time and memory.
    Returns
    -------
    out: scipy.sparse.csr_matrix
        Convolution product of signal by kernel.
    """
    sm, sn = signal.shape

    if type(kernel) is tuple:
        kernel_l, kernel_r = kernel
        km = kernel_l.shape[0]
        kn = kernel_r.shape[1]
        if kernel_l.shape[1] != kernel_r.shape[0]:
            raise ValueError("Kernel factorisation is invalid")
        n_factors = kernel_l.shape[1]
        for f in range(n_factors):
            subkernel_l = sp.diags(
                kernel_l[:, f],
                np.arange(km),
                shape=(sm - km + 1, sm),
                format="dia",
            )
            subkernel_r = sp.diags(
                kernel_r[f, :],
                -np.arange(kn),
                shape=(sn, sn - kn + 1),
                format="dia",
            )
            if f == 0:
                out = (subkernel_l @ signal) @ subkernel_r
            else:
                out += (subkernel_l @ signal) @ subkernel_r
    else:
        km, kn = kernel.shape

        # Sanity checks
        if sp.issparse(kernel):
            raise ValueError("cannot handle kernel in sparse format")
        if not sp.issparse(signal):
            raise ValueError("cannot handle signal in dense format")
        # Check of kernel is constant (uniform)
        constant_kernel = np.nan
        if np.allclose(
            kernel, np.tile(kernel[0, 0], kernel.shape), rtol=1e-08
        ):
            constant_kernel = kernel[0, 0]

        out = sp.csc_matrix((sm - km + 1, sn - kn + 1), dtype=np.float64)

        # Simplified convolution for the special case where kernel is constant:
        if np.isfinite(constant_kernel):
            l_subkernel_sp = sp.diags(
                constant_kernel * np.ones(km),
                np.arange(km),
                shape=(sm - km + 1, sm),
                format="dia",
            )
            r_subkernel_sp = sp.diags(
                np.ones(kn),
                -np.arange(kn),
                shape=(sn, sn - kn + 1),
                format="dia",
            )
            out = (l_subkernel_sp @ signal) @ r_subkernel_sp
        # Convolution code for general case
        else:
            for kj in range(kn):
                subkernel_sp = sp.diags(
                    kernel[:, kj],
                    np.arange(km),
                    shape=(sm - km + 1, sm),
                    format="csr",
                )
                out += subkernel_sp.dot(signal[:, kj : sn - kn + 1 + kj])

    # Set very low pixels to 0
    out.data[np.abs(out.data) < threshold] = 0
    out.eliminate_zeros()

    # Resize matrix to original dimensions
    out = preproc.zero_pad_sparse(
        out, margin_h=(kn - 1) // 2, margin_v=(km - 1) // 2, fmt="csr"
    )
    return out

# Input arguments:
df1= pd.read_table(sys.argv[1],header=0, delimiter="\t")
bank1 = sys.argv[2]
df2= pd.read_table(sys.argv[3],header=0, delimiter="\t")
bank2= sys.argv[4]

print("Number of elements in group1:")
print(len(df1))
print("Number of elements in group2:")
print(len(df2))
print("Sum:")
print(len(df1)+len(df2))

jitter = 3 # total size of jitter around detected coordinates

l1=list(df1['chrom1'])
l2=list(df1['chrom2'])
l3=list(df2['chrom1'])
l4=list(df2['chrom2'])
list_all_chrms=np.unique(l1+l2+l3+l4)

print("Number of chromosomes:")
print(len(list_all_chrms))

ci=0
n_group1=0
n_group2=0
n_group12=0
for c in list_all_chrms :
    ci+=1
    if ci ==1:
        bool_val = True
        mode_val= "w"
    else : 
        bool_val= False
        mode_val= "a"
        
    # 2D coordinates belonging to the same chr
    ka = df1.loc[ ( (df1['chrom1'] == c) & (df1['chrom2'] == c) )]
    kb = df2.loc[ ( (df2['chrom1'] == c) & (df2['chrom2'] == c) )]
    
    # Conversion into sparse objects:
    row_a = np.array(ka['bin1'])
    col_a = np.array(ka['bin2'])
    
    row_b = np.array(kb['bin1'])
    col_b = np.array(kb['bin2'])
    
    mx = max( np.concatenate((row_a, row_b), axis=0) ) 
    my = max( np.concatenate((col_a, col_b), axis=0) ) 
    
    data_a = np.ones( len(row_a) )
    data_b = np.ones( len(row_b) )
    
    csr_a = sp.csr_matrix( (data_a,(row_a,col_a)) ,shape=(mx+jitter,my+jitter))
    csr_b = sp.csr_matrix( (data_b,(row_b,col_b)) ,shape=(mx+jitter,my+jitter))
    
    # adding of jitter with sparse convolution:
    kernel = np.ones( (jitter,jitter) )
    csr_a_j = xcorr2_sparse(csr_a, kernel, threshold=1e-6)
    csr_b_j = xcorr2_sparse(csr_b, kernel, threshold=1e-6)
    
    # Sum sparse object:
    csr_sum = csr_a_j + csr_b_j 
    
    # we assign each coordinate in each of the 3 sub-groups:
    e1=[]    
    e12=[]
    for i in range( len(ka) ):
        xi=row_a[i]
        yi=col_a[i]
        if (csr_sum[xi-int(jitter/2.):xi+int(jitter/2.+1),
                    yi-int(jitter/2.):yi+int(jitter/2.+1)]).sum() <= jitter**2 :
            e1.append(i)
        else : 
            e12.append(i)
        
    ka1 = ka.iloc[e1]
    ka1.to_csv("group1_detected_in_"+bank1+".txt", sep='\t', 
               header=bool_val, index=False, mode=mode_val) 
    n_group1=n_group1+len(ka1)
    
    ka12 = ka.iloc[e12]
    ka12.to_csv("group12_detected_in_"+bank1+"_"+bank2+".txt", sep='\t',
                header=bool_val, index=False, mode=mode_val)
    n_group12=n_group12+len(ka12) 
    
    e2=[]
    for i in range( len(kb) ):
        xi=row_b[i]
        yi=col_b[i]
        if (csr_sum[xi-int(jitter/2.):xi+int(jitter/2.+1),
                    yi-int(jitter/2.):yi+int(jitter/2.+1)]).sum() <= jitter**2 :
            e2.append(i)
        else : 
            e12.append(i)
    kb2 = kb.iloc[e2]
    kb2.to_csv("group2_detected_in_"+bank2+".txt", sep='\t', 
               header=bool_val, index=False, mode=mode_val)
    n_group2=n_group2+len(kb2) 

print("Number of elements only in group1:")    
print(n_group1)  
print("Number of elements only in group2:")
print(n_group2)
print("Number of elements in group12 (intersection):")
print(n_group12)

if len(df1) + len(df2) != n_group1 + n_group2 + 2*n_group12 :
    print("Error: the sum of elements are different.")

# PLot:
# First way to call the 2 group Venn diagram:
venn2(subsets = (n_group1, n_group2, n_group12), 
      set_labels = (bank1, bank2), set_colors=  ('blue', 'pink') )
plt.savefig('venn_diagram_group_loops_'+bank1+'_'+bank2+'.pdf')
#plt.show()
plt.close()
    

    
