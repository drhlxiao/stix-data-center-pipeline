#!/usr/bin/python
"""
    Process ScienceL1 fits file
    used by preview image creation software
    Author: Hualin Xiao (hualin.xiao@fhnw.ch)
    Date: Sep. 1, 2021
"""
import datetime
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from stix.spice import time_utils as sdt
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle



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
        #timebin: detector,pixel, energies

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
    

    def get_time_and_counts(self, emin_sci:int, emax_sci:int, integration_time:float, start=None, end=None):
        """
            find the peak time and integrate counts 
        """

        counts=np.sum(self.spectrogram[:,emin_sci:emax_sci], axis=1)
        if counts.size==0:
            return None, None, None


        if start is None or end is None:
            #use the peak time
            max_idx = np.argmax(counts)
            #max_rate =np.max(counts)/self.timedel[max_idx]
            peak_time_bin=[ self.time_bins_low[max_idx], self.time_bins_high[max_idx]]
            peak_time=(peak_time_bin[0]+peak_time_bin[1])/2.
            #mean time
            if peak_time_bin[1]-peak_time_bin[0] >= integration_time:
                #one time bin requests
                start,end=peak_time_bin
            else:
                start = max(peak_time-integration_time/2, self.time_bins_low[0])
                end = min(peak_time+integration_time/2, self.time_bins_high[-1])
        
        start=max(start, self.time_bins_low[0])
        end=min(end, self.time_bins_high[-1])
        #make sure start time and end valid


        pixel_total_counts=np.sum(self.counts[  (self.time_bins_low >=start) & (self.time_bins_high <= end) ,:,:,emin_sci:emax_sci], axis=(0,3))
        #timebin: detector,pixel, energies

        total_counts=np.sum(pixel_total_counts)
        return start, end, total_counts, pixel_total_counts
    
    def plot_spectrogram(self, ax=None, selection_box=None):
        if not ax:
            _, ax = plt.subplots()
        X, Y = np.meshgrid(self.datetime,
                           np.arange(self.min_ebin, self.max_ebin))
        im = ax.pcolormesh(
                X, Y,
                np.transpose(
                    self.count_rate_spectrogram[:, self.min_ebin:self.max_ebin]
                ))  #pixel summed energy spectrum
        ax.set_yticks(
                self.energies['channel'][self.min_ebin:self.max_ebin:2])
        ax.set_yticklabels(
                self.energy_bin_names[self.min_ebin:self.max_ebin:2])
        fig = plt.gcf()
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('Counts')
        ax.set_title('Count rate spectrogram')
        ax.set_ylabel('Energy range (keV')

        if selection_box:
            erange=selection_box.get('erange',None)
            trange=selection_box.get('trange',None)
            if erange and trange:
                emin_sci, emax_sci=self.energy_to_index(erange[0], erange[1])
                tmin,tmax=sdt.utc2datetime(trange[0]), sdt.utc2datetime(trange[1]) 
                width=mdates.date2num(tmax)-mdates.date2num(tmin)

                rec= Rectangle((tmin, emin_sci), 
                    width, emax_sci-emin_sci)
                    
                ax.add_patch(rec, fill=None,edgecolor='cyan')

                locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
                formatter = mdates.ConciseDateFormatter(locator)
                ax.xaxis.set_major_locator(locator)
                ax.xaxis.set_major_formatter(formatter)




        #ax.set_xlabel(f"Time")
        return ax


    def get_time_ranges_for_imaging(self, flare_time_ranges,  imaging_energies,  min_counts=3000, integration_time = 60):
        """
        determine time ranges for imaging
            elow_keV: float
                energy range lower limit
            ehigh_keV: float
                energy range upper limit
            min_counts: int
                minimum counts per bin
            """
        boxes=[]
        sci_energy_ranges=[]

        for energy_range in imaging_energies:
            elow_sci, emax_sci=self.energy_to_index(energy_range[0], energy_range[1])
            if elow_sci is None or emax_sci is None:
                continue
            sci_energy_ranges.append([elow_sci, emax_sci])

        time_ranges = []

        for i, sci_range in enumerate(sci_energy_ranges):

            start, end, total_counts, pixel_counts =self.get_time_and_counts( sci_range[0], sci_range[1], integration_time)
            if start is None:
                continue
            if total_counts>min_counts:
                boxes.append({
                    'total_counts':  total_counts,
                    'pixel_counts':  pixel_counts,
                    'energy_range_sci': sci_range,
                    'energy_range_keV': imaging_energies[i],
                    'unix_time_range': [start,  end],
                    'utc_range': [sdt.unix2utc(start),  sdt.unix2utc(end)]})
            else:
                #this might be a small flare, take the whole time period
                for flare in flare_time_ranges:
                    #create images for flaring time
                    if flare['peak']:
                        #try again
                        start, end, total_counts, pixel_counts =self.get_time_and_counts(sci_range[0], sci_range[1], integration_time, flare['start'],flare['end'])
                        if total_counts>min_counts:
                            boxes.append({
                                'total_counts':  total_counts,
                                'pixel_counts':  pixel_counts,
                                'energy_range_sci': sci_range,
                                'energy_range_keV': imaging_energies[i],
                                'unix_time_range': [start,  end],
                                'utc_range': [sdt.unix2utc(start),  sdt.unix2utc(end)]})




        return boxes

                




