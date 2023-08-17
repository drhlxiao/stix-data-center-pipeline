#!/usr/bin/python
import numpy as np
'''
   data models
   spectrogram and spectrum

    - dead time correction
    - transmission correction
    - energy binning correction
    Author: Hualin Xiao
    Date: Aug. 29, 2021

'''

DETECTOR_GROUPS = [[1, 2], [6, 7], [5, 11], [12, 13], [14, 15], [10, 16],
                   [8, 9], [3, 4], [31, 32], [26, 27], [22, 28], [20, 21],
                   [18, 19], [17, 23], [24, 25], [29, 30]]


class Spectrogram(object):
    '''
       Spectrogram formatter
    '''

    def __init__(self):
        self.counts = []
        self.time_bins = []
        self.ebins = []

    def fill(self, ebin_low, ebin_up, time, counts):
        if time not in self.time_bins:
            self.time_bins.append(time)
            self.counts.append([0] * 32)
        itbin = self.time_bins.index(time)
        ebin = (ebin_low, ebin_up)
        if ebin not in self.ebins:
            self.ebins.append(ebin)
        iebin = self.ebins.index(ebin)

        self.counts[itbin][iebin] += counts

    def get_spectrogram(self, dtype='numpy.array'):
        spec = np.array(self.counts)
        num_bins = len(self.ebins)
        sp = spec[:, 0:num_bins]  # max 32 bins, slice the array
        spectrogram = sp.tolist() if dtype == 'list' else sp
        return self.time_bins, self.ebins, spectrogram


class EnergySpectrum(object):
    '''
        do corrections for x-ray flux spectrum

    '''

    def __init__(self):
        self.spectrum = []
        self.duration = []
        self.start_time = 0
        self.end_time = 0
        self.tau = 0.4e-6

    def fill(self, unix_time, ebin_low, ebin_up, bin_counts, triggers):

        detector_id = pixel_id // 12

    def fill_one_detector(self, unix_time, time_bin, detector_id, ebin_low,
                          ebin_up, bin_counts, triggers):

        pass

    def __sub__(self, B):
        pass
