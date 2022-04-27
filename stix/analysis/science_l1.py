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

        if self.is_time_bin_shifted(self.T0_unix) and len(self.timedel)>1:
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
    def energy_to_index(self, elow_keV, ehigh_keV):
        ebin_min=None
        ebin_max=None
        for  e in self.energies:
            #print(e)
            if e[1]<=elow_keV:
                ebin_min=e[0]
            if e[2]<=ehigh_keV:
                ebin_max=e[0]
        if ebin_min is None or ebin_max is None:
            return None,None
        ebin_max+=1

        return ebin_min, ebin_max

    def make_spectra(self, pixel_counts=None):
        #print('spectrogram...')
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

        self.time_bins_low=self.time-0.5*self.timedel+self.T0_unix
        self.time_bins_high=self.time+0.5*self.timedel+self.T0_unix

    def get_total_counts(self, emin_sci:int, emax_sci:int, unix_start, unix_end):
        counts=np.sum(self.spectrogram[:,emin_sci:emax_sci], axis=1)
        if unix_end is None or unix_end is None:
            return np.sum(counts)
        total_counts=np.sum(counts[ (self.time_bins_low >=unix_start) & (self.time_bins_high <= unix_end)])
        return total_counts


    def get_time_ranges_for_imaging(self, imaging_energies, flare_unix_time_ranges, min_counts=3000, integration_time = 60, time_step=300):
        """
        determine time ranges for imaging
            flare_unix_time_ranges: list
                flare time ranges
            elow_keV: float
                energy range lower limit
            ehigh_keV: float
                energy range upper limit
            min_counts: int
                minimum counts per bin
            """


        boxes=[]
        if not flare_unix_time_ranges:
            return []

        sci_energy_ranges=[]
        for energy_range in imaging_energies:
            elow_sci, emax_sci=self.energy_to_index(energy_range[0], energy_range[1])
            if elow_sci is None or emax_sci is None:
                sci_energy_ranges.append(None)
            sci_energy_ranges.append([elow_sci, emax_sci])

        time_ranges = []


        duration = self.time[-1]+ 0.5* self.timedel[-1]-(self.time[0]+ 0.5* self.timedel[0])
        if duration <= integration_time or len(self.timedel)==1:
            # duration too short or only one time bin
            start = self.time[0]- 0.5* self.timedel[0] + self.T0_unix
            end = self.time[-1]+ 0.5* self.timedel[-1] +self.T0_unix
            time_ranges= [[start, end]]
        else:
            for flare_time in flare_unix_time_ranges:
                #only select flaring times
                start, peak, end = flare_time
                time_ranges.append( [start, end]) #flare  whole time range

                peak= peak if peak is not None else start
                if peak is not None:
                    #every time_step
                    t=peak
                    while t-integration_time/2 >= start :
                        time_ranges.append( [t-integration_time/2., t+integration_time/2.])
                        t = t - time_step
                t=peak
                print('Peak time:', sdt.unix2utc(peak))
                while t+integration_time/2. <= end :
                    time_ranges.append( [t-integration_time/2., t+integration_time/2.])
                    t = t + time_step
        

        for flare_time in time_ranges:
            #only select flaring times
            start, end = flare_time

            total_counts = [0]*len(imaging_energies)
            
            for i, sci_range in enumerate(sci_energy_ranges):
                if sci_range:
                    total_counts[i]= self.get_total_counts(sci_range[0], sci_range[1], start, end)
                    #don't make images if count rate too low 
                #both energies don't have counts
            boxes.append({
                'total_counts':  total_counts,
                'counts_enough':  [bool(x>min_counts) for x in total_counts],
                    'energy_range_sci': sci_energy_ranges,
                    'energy_range_keV': imaging_energies,
                    'unix_time_range': [start,  end],
                    'utc_range': [sdt.unix2utc(start),  sdt.unix2utc(end)]})
        return boxes

                




