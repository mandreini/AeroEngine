# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:34:17 2019

@author: Matt
"""
import math
import matplotlib.pyplot as plt
import numpy as np
import engine
import characteristics
#import sankeyview

def get_log(svals, Tvals, lind):
    # This takes 2 (x, y) coordinates (in this case, T and s) and gets a 
    # logarithmic fit between them to generate a small curve. This is just
    # a hack to get a curve for the T-s diagram.
    isenfit = np.polyfit(svals[lind:lind+2], np.log(Tvals[lind:lind+2]), 1)
    xs = np.linspace(svals[lind], svals[lind+1])
    ys = np.exp(isenfit[1]) * np.exp(isenfit[0] * xs)
    xs, ys = list(xs), list(ys)        
    
    svals = svals[:lind] + xs + svals[lind+1:]
    Tvals = Tvals[:lind] + ys + Tvals[lind+1:]
    
    return svals, Tvals 
    
def plot_all(dryeng, abeng, itbeng):
    # plot and label the T-s diagram for all 3 configurations for AET3
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #dry config  
    if dryeng is not None:
        dry_Tvals = dryeng[0]
        dry_ds = dryeng[1]
        dry_svals = [sum(dry_ds[:i]) for i in range(len(dry_ds))]
#        dry_coords = get_log(dry_svals, dry_Tvals, 3)
        dry_coords = (dry_svals, dry_Tvals)
        plt.plot(dry_coords[0], dry_coords[1], label=dryeng[2], color='r')
        
        ax.annotate('0 [%i, %i]' % (round(dry_Tvals[0]), round(dry_svals[0])), xy=(dry_svals[0]+25, dry_Tvals[0]-50), textcoords='data')
        ax.annotate('2 [%i, %i]' % (round(dry_Tvals[1]), round(dry_svals[1])), xy=(dry_svals[1]+35, dry_Tvals[1]), textcoords='data')
        ax.annotate('25 [%i, %i]' % (round(dry_Tvals[2]), round(dry_svals[2])), xy=(dry_svals[2]+35, dry_Tvals[2]-15), textcoords='data')
        ax.annotate('3 [%i, %i]' % (round(dry_Tvals[3]), round(dry_svals[3])), xy=(dry_svals[3]+75, dry_Tvals[3]-15), textcoords='data')
        ax.annotate('4 [%i, %i]' % (round(dry_Tvals[4]), round(dry_svals[4])), xy=(dry_svals[4], dry_Tvals[4]+50), textcoords='data')
        ax.annotate('            45 \n [%i, %i]' % (round(dry_Tvals[5]), round(dry_svals[5])), xy=(dry_svals[5]-150, dry_Tvals[5]-25), textcoords='data')
        ax.annotate('5 [%i, %i]' % (round(dry_Tvals[7]), round(dry_svals[7])), xy=(dry_svals[7]+35, dry_Tvals[7]-25), textcoords='data')
        ax.annotate('8 [%i, %i]' % (round(dry_Tvals[9]), round(dry_svals[9])), xy=(dry_svals[9]+25, dry_Tvals[9]), textcoords='data')
        
    
    if abeng is not None:
        ab_Tvals = abeng[0]
        ab_ds = abeng[1]
        ab_svals = [sum(ab_ds[:i]) for i in range(len(ab_ds))]
#        ab_coords = get_log(ab_svals, ab_Tvals, 8)
#        ab_coords = get_log(ab_coords[0], ab_coords[1], 3) # not working?
        ab_coords = (ab_svals, ab_Tvals)
        plt.plot(ab_coords[0], ab_coords[1], label=abeng[2], color='b')
        
        ax.annotate('7 (ab) \n[%i, %i]' % (np.round(ab_Tvals[8]), np.round(ab_svals[8])), xy=(ab_svals[8]+15, ab_Tvals[8]), textcoords='data')
        ax.annotate('8 (ab) \n[%i, %i]' % (np.round(ab_Tvals[9]), np.round(ab_svals[9])), xy=(ab_svals[9]+15, ab_Tvals[9]), textcoords='data')
#        
    if itbeng is not None:
        itb_Tvals = itbeng[0]
        itb_ds = itbeng[1]
        itb_svals = [sum(itb_ds[:i]) for i in range(len(itb_ds))]
#        itb_coords = get_log(itb_svals, itb_Tvals, 7)
#        itb_coords = get_log(itb_svals, itb_Tvals, 3)
        itb_coords = (itb_svals, itb_Tvals)
        plt.plot(itb_coords[0], itb_coords[1], label=itbeng[2], color='g')
        
        ax.annotate('47 (itb) [%i, %i]' % (np.round(itb_Tvals[6]), np.round(itb_svals[6])), xy=(itb_svals[6]-200, itb_Tvals[6]+25), textcoords='data')
        ax.annotate('     5 (itb) \n[%i, %i]' % (np.round(itb_Tvals[7]), np.round(itb_svals[7])), xy=(itb_svals[7]-215, itb_Tvals[7]-50), textcoords='data')
        ax.annotate('8 (itb) \n[%i, %i]' % (np.round(itb_Tvals[9]), np.round(itb_svals[9])), xy=(itb_svals[9]-80, itb_Tvals[9]-50), textcoords='data')

    mdry = (dry_Tvals[4]-dry_svals[3])/(dry_svals[4]-dry_svals[3])
    mab = (ab_Tvals[8]-ab_svals[7])/(ab_svals[8]-ab_svals[7])
    mitb = (itb_Tvals[6]-itb_Tvals[5])/(itb_svals[6]-itb_svals[5])

    plt.legend(loc='upper left')
    plt.xlabel('s [J/K]')
    plt.ylabel('T [K]')
    plt.title('Comparison of T-s diagram of all three engines')
    
    return mdry, mab, mitb

def generate_sankey(effs, fname):
    pass
