#!/usr/bin/python
"""
    Process ScienceL1 fits file
    Author: Hualin Xiao (hualin.xiao@fhnw.ch)
    Date: Sep. 1, 2021
"""
import datetime
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from stix.spice import time_utils as sdt



class ScienceData():
    """
      Retrieve science data from stix data center or load fits file from local storage

    """
    def __init__(self, request_id=None, fname=None):
        self.fname = fname
        self.data_type=None
        if not fname:
            raise Exception("FITS filename not specified")
        self.request_id = request_id
        self.time_shift_applied=0
        self.hdul = fits.open(fname)
        self.energies=[]
        #self.read_data()

    def read_fits(self, light_time_correction=True):
        """
            Read data  L1 compression level  FITS files
            Parameters
            ---------------------
            light_time_correction: boolean
                Correct light time difference
        """

        self.data = self.hdul['DATA'].data
        self.T0_utc = self.hdul['PRIMARY'].header['DATE_BEG']
        self.counts= self.data['counts']

        self.light_time_del= self.hdul['PRIMARY'].header['EAR_TDEL']
        self.light_time_corrected=light_time_correction

        self.T0_unix = sdt.utc2unix(self.T0_utc)
        self.triggers = self.data['triggers']
        self.rcr = self.data['rcr']

        self.timedel = self.data['timedel']
        self.time = self.data['time']

        if self.is_time_bin_shifted(self.T0_unix):
            self.timedel = self.timedel[:-1]
            self.time = self.time[1:]
            print('Shifted time bins have been corrected automatically!')
            if self.data_type=='ScienceL1':
                self.counts= self.counts[1:, :, :, :]
                self.triggers = self.triggers[1:, :]
                self.rcr = self.rcr[1:]
            elif self.data_type=='Spectrogram':
                self.counts= self.counts[1:, :]
                self.triggers = self.triggers[1:]
                #self.rcr = self.rcr[1:]

        self.request_id = self.hdul['CONTROL'].data['request_id']
        
        self.time_shift_applied=0 if light_time_correction else self.light_time_del
        self.datetime = [
            sdt.unix2datetime(self.T0_unix + x + y * 0.5 + self.time_shift_applied)
            for x, y in zip(self.time, self.timedel)
        ]

        self.duration = self.time[-1] - self.time[0] + (self.timedel[0] +
                                                        self.timedel[-1]) / 2

        self.energies = self.hdul['ENERGIES'].data

        self.energy_bin_names = [
            f'{a} - {b}'
            for a, b in zip(self.energies['e_low'], self.energies['e_high'])
        ]
        self.energy_bin_mask = self.hdul["CONTROL"].data["energy_bin_mask"]

        ebin_nz_idx =self.energy_bin_mask.nonzero()
        self.max_ebin = np.max(ebin_nz_idx)  #indices of the non-zero elements
        self.min_ebin = np.min(ebin_nz_idx)

        self.ebins_mid = [
            (a + b) / 2.
            for a, b in zip(self.energies['e_low'], self.energies['e_high'])
        ]
        self.ebins_low, self.ebins_high = self.energies[
            'e_low'], self.energies['e_high']

        if self.data_type=='ScienceL1':
            self.pixel_counts=self.counts
            self.pixel_count_rates= self.pixel_counts/self.timedel[:,None,None, None]
            self.trigger_rates = self.triggers / self.timedel[:, None] 
        elif self.data_type=='Spectrogram':
            self.count_rates= self.counts/self.timedel[:,None]
            self.trigger_rates = self.triggers / self.timedel


    def is_time_bin_shifted(self, unix_time):
        """
            Time bins are shifted in the data collected before 2021-12-09 due a bug in the flight software

            Check if time bin is shifted in L1 data
            Parameters
                unix_time: float 
            Returns
                is_shifted: bool
                    True if time bin is shifted else False
        """

        return (unix_time < sdt.utc2unix('2021-12-09T14:00:00'))



    @classmethod
    def from_fits(cls, filename):
        """
        factory class
        Arguments
        filename: str
            FITS filename
        """
        request_id = None
        return cls(request_id, filename)

    def get_energy_range_slicer(self, elow, ehigh):
        sel=[]
        i=0
        for a, b in zip(self.energies['e_low'], self.energies['e_high']):
            if a>=elow  and b<=ehigh:
                sel.append(i)
            i+=1
        return slice(min(sel),max(sel))



    def __getattr__(self, name):
        if name == 'data':
            return self.hdul
        elif name == 'type':
            return self.hdul.get('data_type', 'INVALID_TYPE')
        elif name == 'filename':
            return self.fname

    def get_data(self):
        return self.hdul


class ScienceL1(ScienceData):
    """
    Tools to analyze L1 science data
    """

    def __init__(self, reqeust_id,fname):
        super().__init__(reqeust_id, fname)
        self.data_type='ScienceL1'
        self.pixel_count_rates=None
        self.correct_pixel_count_rates=None
        self.read_fits()
        self.make_spectra()

    def make_spectra(self, pixel_counts=None):
        print('spectrogram...')
        if pixel_counts is None:
            pixel_counts=self.pixel_counts
        self.spectrogram = np.sum(pixel_counts, axis=(1, 2))
        self.count_rate_spectrogram = self.spectrogram / self.timedel[:, np.
                                                                      newaxis]
        self.spectrum = np.sum(pixel_counts, axis=(0, 1, 2))
        self.mean_pixel_rate_spectra = np.sum(self.pixel_counts,
                                              axis=0) / self.duration
        self.mean_pixel_rate_spectra_err = np.sqrt(
            self.mean_pixel_rate_spectra) / np.sqrt(self.duration)
        #sum over all time bins and then divide them by the duration, counts per second
    def get_time_bins_for_imaging(self, elow_keV=6, ehigh_keV=10, min_counts=5000, min_duration = 60 , signal_unix_time_range=[-np.inf,np.inf]):
        """
         automatic select time ranges for imaging
         Arguments
         ----
         signal_unix_time_range: list
            signal unix time range
         min_counts: int
            minimal counts per bin
        min_duration: int
            minimal time bin 
        """
        time_ranges=[]
        num_tbins=len(self.time)
        last_tbin = [self.time[0]- 0.5* self.timedel[0], self.time[0]+ 0.5* self.timedel[1]]
        
        if num_tbins==1:
            return np.array(last_tbin), 0
        ebin_min=None
        ebin_max=None
        for  e in self.energies:
            #print(e)
            if e[1]<=elow_keV:
                ebin_min=e[0]
            if e[2]<=ehigh_keV:
                ebin_max=e[0]

        if ebin_min is None or ebin_max is None:
            return None,0
        #print('Science bins:',ebin_min, ebin_max)
        ebin_max+=1
        #print(self.spectrogram.shape)
        counts=np.sum(self.spectrogram[:,ebin_min:ebin_max], axis=1)
        cnt_sum=0
        
        
        
        cnt_sum=0
        summed_cnts=[]
        #print(counts.shape)
        for i, c in enumerate(counts):
            cnt_sum+=c
            if self.time[i] - 0.5*self.timedel[i] < signal_unix_time_range[0]:
                continue

            this_tbin=[last_tbin[1], self.time[i]+0.5*self.timedel[i]]
            if cnt_sum>=min_counts and this_tbin[1]-this_tbin[0] >= min_duration:
                #t1=self.time[i]
                time_ranges.append(this_tbin)
                summed_cnts.append(cnt_sum)
                last_tbin= this_tbin
            if self.time[i] + 0.5*self.timedel[i] > signal_unix_time_range[1]:
                break
                
        return np.array(time_ranges)+self.T0_unix, summed_cnts


