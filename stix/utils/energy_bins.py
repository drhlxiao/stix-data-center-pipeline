#!/usr/bin/python
'''
   Any routines related to energy bin conversions are placed here
   Author: Hualin Xiao
   Date: Aug 26, 2021
'''
EBINS = [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150, 'inf']
def to_keV(low_bin:int, up_bin:int):
    if low_bin<0 or low_bin>=32 or up_bin>=32 or up_bin<0:
        return 'Invalid Ebin'
    elif up_bin == 31:
        return f'{EBINS[low_bin]}  keV - inf';
    return f'{EBINS[low_bin]}  -  {EBINS[up_bin + 1]}  keV';

