#!/usr/bin/python
'''
   Any routines related to energy bin conversions are placed here
   Author: Hualin Xiao
   Date: Aug 26, 2021
'''
import numpy as np
EBINS = [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150, np.inf]
ELUT_ENERGY_BINS= [
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36,
    40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150
]
EBINS_LOW=np.array(EBINS)[:-1]
EBINS_HIGH=np.array(EBINS)[1:]

def to_keV(low_bin:int, up_bin:int):
    if low_bin<0 or low_bin>=32 or up_bin>=32 or up_bin<0:
        return 'Invalid Ebin'
    elif up_bin == 31:
        return f'{EBINS[low_bin]}  keV - inf';
    return f'{EBINS[low_bin]}  -  {EBINS[up_bin + 1]}  keV';

def get_sci_energy_bins():
    return EBINS

def keV2sci(elow_keV, ehigh_keV):
    sel_bins=[i for i in range(32)   if EBINS_LOW[i]>=elow_keV and EBINS_HIGH>=ehigh_keV]
    try:
        min(sel_bins), max(sel_bins)
    except ValueError:
        return None, None
def sci2keV(elow_sci, ehigh_sci):
    return EBINS_LOW[elow_sci], EBINS_HIGH[ehigh_sci]


def get_emask_energy_bins(emask):
    ebins=[]
    for i in range(32):
        if emask & (1 << i) != 0:
            ebins.append(i);
    names=[]
    sci_edges=[]
    for i in range(len(ebins) - 1):
        begin = ebins[i]
        end = ebins[i + 1]
        sci_edges.append([begin,end])
        if end == 32:
            names.append(f'{EBINS[begin]} keV –⁠ inf')
        elif end < 32:
            names.append(f'{EBINS[begin]}  – {EBINS[end]} keV')
        else:
            names.append('')
    return {'names':names, 'sci_bin_edges':sci_edges}



def get_corrected_energy_bins(utc, sci_elow, sci_eup):
    #use real elut to correct energy bins
    return (EBINS[sci_elow], EBINS[sci_eup+1])
def get_elut_energy_bins():
    return ELUT_ENERGY_BINS
