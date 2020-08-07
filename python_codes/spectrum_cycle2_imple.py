#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 18:21:09 2020
@author: axel KournaK 
Computute directly on the scatter plot (<1Mb) the loess and plot 
"""
import numpy as np
import scipy
import time
import itertools
import matplotlib.pylab as plt
import json
import pandas as pd
import matplotlib.gridspec as gridspec
import sys 
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
from skmisc.loess import loess

times_points = ["0h","2h15","2h30","3h","4h","8h"] 
times_points.reverse()

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

# another coding of colours 
from colour import Color
pink = Color("pink")
royalblue = Color("royalblue")
colors = list(royalblue.range_to(Color("pink"), len(times_points) +1 ))

# Plot of several time points: 
step=10000
list_all_sizes = range(0,1000000,step) 
i=0

plt.axvline(x=130000,linestyle='dotted',color='black', linewidth=0.5)
plt.axvline(x=260000,linestyle='dotted',color='black', linewidth=0.5)
plt.axvline(x=390000,linestyle='dotted',color='black', linewidth=0.5)
plt.axvline(x=455000,linestyle='dotted',color='black', linewidth=0.5)

for T in times_points:
    print(T)
    i+=1
    if T=="0h":
        txt_file = "GSM3909703_TB-HiC-Dpn-R2-T0_hg19.1000.multires.cool.res8000.cool.res.1Mb/loops_quant.txt"
        
    if T=="2h":
        txt_file = "GSM3909698_TB-HiC-Dpn-R2-T2_hg19.1000.multires.cool.res8000.cool.res.1Mb/loops_quant.txt"
    
    if T=="2h15":
        txt_file = "GSM3909697_TB-HiC-Dpn-R2-T225_hg19.1000.multires.cool.res8000.cool.res.1Mb/loops_quant.txt"
    
    if T=="2h30":
        txt_file = "GSM3909696_TB-HiC-Dpn-R2-T25_hg19.1000.multires.cool.res8000.cool.res.1Mb/loops_quant.txt"
    
    if T=="3h":
        txt_file = "GSM3909694_TB-HiC-Dpn-R2-T3_hg19.1000.multires.cool.res8000.cool.res.1Mb/loops_quant.txt"
    
    if T=="4h":
        txt_file = "GSM3909691_TB-HiC-Dpn-R2-T4_hg19.1000.multires.cool.res8000.cool.res.1Mb/loops_quant.txt"
    
    if T=="8h":
        txt_file = "GSM3909686_TB-HiC-Dpn-R2-T8_hg19.1000.multires.cool.res8000.cool.res.1Mb/loops_quant.txt"
    
    if T=="10h":
        txt_file = "GSM3909684_TB-HiC-Dpn-R2-T10_hg19.1000.multires.cool.res8000.cool.res.1Mb/loops_quant.txt"
    
    if T=="12h":
        txt_file = "GSM3909682_TB-HiC-Dpn-R2-T12_hg19.1000.multires.cool.res8000.cool.res.1Mb/loops_quant.txt"
    
    print(i,2.0/i)
    color1=lighten_color('royalblue', i*0.2)
    color1= colors[i].hex_l
#    color2=  color=lighten_color(color1, (i-0.5)*0.2)
    #plt.plot(list_all_sizes, list_all_scores, color=color1, label="T0") 
#    plt.plot(list_all_sizes, list_all_scores,'o', color = color1, label=T)
    
    #-------------------
    threshold_size = 1000000 # maximal size for the computation 
    df = pd.read_csv(txt_file,header=0, delimiter="\t")   # txt file
    df = df[(abs(df['start2']-df['start1'])<threshold_size) ]
    print(len(df))
    
    list_all_sizes = abs (df['start2']-df['start1'] )
    list_all_scores = df['score']
    
    list_all_sizes = list_all_sizes[~np.isnan(list_all_scores)]
    list_all_scores = list_all_scores[~np.isnan(list_all_scores)]
    
    list_all_sizes = np.array(list_all_sizes)
    list_all_scores = np.array(list_all_scores)
    
    list_all_scores = list_all_scores[np.argsort(list_all_sizes)]
    list_all_sizes.sort()
    
    #-------------------
    x = list_all_sizes
    y = list_all_scores
    #plt.plot(x,y,'o')
    
    l = loess(x,y, span=0.15)   # do not take a big percentage
    l.fit()
    pred = l.predict(x, stderror=True)  #  time consuming step!
    conf = pred.confidence()
    
    lowess = pred.values
    ll = conf.lower
    ul = conf.upper
    
    #------------- plot: 
    plt.plot(x, lowess, color = color1, label=T)
    plt.fill_between(x,ll,ul,alpha=.33, color = color1)
    plt.show()
    
    plt.xlabel("Size between two cohesin sites")
    plt.ylabel("Median Loop score")
    
    plt.xlim([0,threshold_size])
    plt.legend()
    plt.title("Loop Spectrum from mitotic to interphase")
    plt.grid()

        
plt.savefig('Loop_score_cohesin_human_types_ALL_span01.pdf')
