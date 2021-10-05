#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find photo-peaks in energy spectra and determine gains and offsets using linear fit
Author: Hualin Xiao (hualin.xiao@fhnw.ch)
Date: Sept. 23, 2021

"""

import os
import sys
import numpy as np
import time
from array import array
from datetime import datetime
from stix.core import stix_datatypes as sdt
from stix.core import mongo_db as db
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from pathlib import Path
from matplotlib import pyplot as plt
from lmfit import Model
from lmfit.models import GaussianModel, LinearModel
from scipy.optimize import curve_fit

from matplotlib.backends.backend_pdf import PdfPages
#SPECTRUM_MODEL_DATA_FILE=Path(__file__).parent.parent/ 'data' / 'ExpSpecModel.npz'

PHOTO_PEAKS_POS = np.array([30.85, 35.13, 81])

mdb = db.MongoDB()

def interp(xvals, yvals, xnew):
    # x y define orignal points
    # xnew interpolated data points
    f2 = interp1d(xvals, yvals, bounds_error=False, fill_value=0)
    y = f2(xnew)
    return y

class PeakFitting(object):

    @staticmethod
    def constant(x,const):
        return const
    @staticmethod
    def linear(x,a, b):
        return a*x+b

    @staticmethod
    def gaussian(x,  amp, cen, wid):
    #    """1-d gaussian: gaussian(x, amp, cen, wid)"""
        return amp*np.exp(-(x-cen)**2 / (2*wid**2))
    @staticmethod
    def asym_gaussian(x, amp, cen, wid, B, C):
        psup = (x >= cen).nonzero()
        pinf = (x < cen).nonzero()
        T_i = np.arange(len(x))
        T_i= (np.exp(B*(x-cen-C)))*(1-np.exp(-(x-cen - C)**2 / (2*wid**2)))
        T_i[psup] = 0.
        return amp*np.exp(-(x-cen)**2 / (2*wid**2)) + amp*T_i
    @staticmethod
    def fit_asym_gaus(x,y, offset=0.1, bw=1):
        moda = Model(PeakFitting.asym_gaussian)+Model(PeakFitting.constant)
        max_counts = np.max(y)  # 80 keV max counts
        B = 0.55
        C = 0.5
        cen = x[y == max_counts][0]
        pars = moda.make_params()
        pars.add('B', value=B, min=0, max=B + 5)
        pars.add('C', value=C, min=0, max=C + 5)
        pars.add('amp', value=1, min=0, max=max_counts)
        pars.add('wid', value=2, min=0, max=5)
        pars.add('cen', value=cen, min=cen-5*bw, max=cen+2*bw)
        pars.add('const', value=offset, min=0, max=0.2*np.max(y))
        result = moda.fit(y, pars, x=x)
        return result

    @staticmethod
    def fit_gaus(x,y, const, amp, cen, wid, bw):
        moda = Model(PeakFitting.gaussian) + Model(PeakFitting.constant)
        pars = moda.make_params()
        pars.add('amp', value=amp, min=0, max=np.max(y))
        pars.add('cen', value=cen, min=cen - 2*bw, max=cen + 2*bw)
        pars.add('wid', value=wid, min=0.1, max=3)
        pars.add('const', value=const, min=0., max=1.1*np.median(y))
        result = moda.fit(y, pars, x=x)
        return result


    @staticmethod
    def fit_double_gaus(x, y,  const,  a_amp, a_cen, a_wid, b_amp, b_cen, b_wid, bw):
        #0, 1, peak1_x, 2, 1, peak2_pos, 1
        moda = Model(PeakFitting.gaussian, prefix='a_') + Model(PeakFitting.gaussian, prefix='b_')+Model(PeakFitting.constant)
        pars = moda.make_params()
        print(a_cen,b_cen)
        pars.add('a_cen', value=a_cen, min=a_cen-2*bw,max=a_cen+2*bw)
        pars.add('b_cen', value=b_cen, min=b_cen - 2*bw, max=b_cen + 2*bw)
        pars.add('a_amp', value=a_amp, min=0, max=np.max(y))
        pars.add('b_amp', value=b_amp, min=0, max=np.max(y))
        pars.add('a_wid', value=a_wid, min=0.1, max=3)
        pars.add('b_wid', value=b_wid, min=0.1, max=3)
        pars.add('const', value=const, min=0., max=1.1*np.median(y))
        result = moda.fit(y, pars, x=x)
        return result
    @staticmethod
    def fit_linear(x, y,yerr):
        moda =  LinearModel()
        #pars = moda.make_params()
        try:
            mean_yerr=np.mean(yerr)
            result = moda.fit(y, x=x, weights=1./yerr)
        except:
            print(x,y, yerr)
            raise
        return result


def select_near_max(x, y, xmin,xmax, left, right):
    """
     get spectrum near the maximum  in the given x range
     Parameters

     x: np.ndarray
     y: np.ndarray
     xmin: float
     xmax: float
        range 
     
     left: float
     right: float
        left and right margin
    """
    _x=np.copy(x)
    _y=np.copy(y)
    _y[(x<xmin)|(x>xmax)]=0
    max_y=np.max(_y)
    max_pos=_x[_y==max_y][0]
    cut= (x > max_pos -left) & (x<max_pos +right)
    return x[cut],y[cut], max_pos

    
class Calibration(object):
    """
    Fit calibration spectrum using ECC methoMEAN_ENERGY_CONVERSION_FACTORd
    """
    FIT_RANGE = (252, 448)
    MEAN_ENERGY_CONVERSION_FACTOR=2.31
    SLOPS_LIMITS = (2.15, 2.45)
    OFFSETS_LIMITS = (205, 245)
    MIN_COUNTS = 100
    DEFAULT_OUTPUT_DIR = '/data/calibration/'
    STEPS=10000
    FIT_RANGE_SMALL_MARGIN = 5
    FIT_RANGE_BIG_MARGIN= 15
    PEAK_MIN_COUNTS = 50
    ELUT_ENERGIES = [
        4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36,
        40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150
    ]
    
    def __init__(self):
        #self.pixel_fspec_models, self.default_slopes, self.default_offsets = \
        #    self.load_spectrum_models(SPECTRUM_MODEL_DATA_FILE)
        pass


    #x=np.linspace(15,90,500)
    #plt.plot(x,fspec_model(x))
    #plt.show()

    def compute_elut(self,offset, slope):
        elut = []
        for det in range(0, 32):
            for pix in range(0, 12):
                p0 = offset[det][pix]
                p1 = slope[det][pix]
                # print(det, pix, p0, p1)

                if p0 > 0 and p1 > 0:
                    row = [det, pix, p0, p1]
                    Elows = [int(4 * (p0 + p1 * x)) for x in Calibration.ELUT_ENERGIES]
                    row.extend(Elows)
                    elut.append(row)
        return elut



    def get_subspec(self, spectrum, start_ch, bin_width, only_fit_range=True):
        """
        reconstruct subspectrum
        Args
        spectrum: List of counts
        start_ch: start ADC channel
        bin_width:  number of summed ADC channels per bin
        xmin, xmax:  limits to slice the spectrum
        Returns
        x: a numpy array indicates ADC channels
        y: a numpy array indicates counts of spectrum
        """
        spec_len = len(spectrum)
        spectrum = np.array(spectrum)
        _x = np.arange(spec_len) * bin_width + 0.5 * bin_width + start_ch
        if only_fit_range:
            clip = (_x > Calibration.FIT_RANGE[0]) & (_x < Calibration.FIT_RANGE[1])
            return _x[clip], spectrum[clip]
        return _x, spectrum


    def find_solution(self, detector, pixel, sbspec_id,spectrum, start_ch, bin_width, pdf=None):
        """
        Fit spectrum using a spectrum model
        Arguments
        spec: spectrum counts
        start_ch: start channel
        bin_width: number of summed ADC bins
        model:
          An 2d array indicate x, y of the theoretical spectrum. Use ECC model if it is None
        """
        spec_x, spec_y = self.get_subspec(spectrum, start_ch, bin_width) #real spectra
        peak1_max_counts = np.max(spec_y)

        peak1_pos = spec_x[spec_y==peak1_max_counts][0]
        #peak 1 has the maximum counts always
        # find the peak with highest counts in the predefined range
        peak1_2_mean_shift=(PHOTO_PEAKS_POS[1]-PHOTO_PEAKS_POS[0]) * Calibration.MEAN_ENERGY_CONVERSION_FACTOR
        peak1_3_mean_shift=(PHOTO_PEAKS_POS[2]-PHOTO_PEAKS_POS[0]) * Calibration.MEAN_ENERGY_CONVERSION_FACTOR

        peak1_x, peak1_y, peak1_pos=select_near_max(spec_x,spec_y, 
                peak1_pos-Calibration.FIT_RANGE_SMALL_MARGIN, peak1_pos+Calibration.FIT_RANGE_SMALL_MARGIN, 
                Calibration.FIT_RANGE_SMALL_MARGIN, Calibration.FIT_RANGE_SMALL_MARGIN)
                
        peak2_x, peak2_y, peak2_pos=select_near_max(spec_x,spec_y,
                peak1_pos + 0.5*peak1_2_mean_shift  , peak1_pos+1.3*peak1_2_mean_shift,  
                Calibration.FIT_RANGE_SMALL_MARGIN, Calibration.FIT_RANGE_SMALL_MARGIN)

        peak12_x, peak12_y, peak12_pos=select_near_max(spec_x,spec_y,
                peak1_pos-Calibration.FIT_RANGE_SMALL_MARGIN, peak2_pos+Calibration.FIT_RANGE_BIG_MARGIN, 
                Calibration.FIT_RANGE_SMALL_MARGIN, peak2_pos-peak1_pos+Calibration.FIT_RANGE_BIG_MARGIN)


        low_e_result=PeakFitting.fit_double_gaus(peak1_x,peak1_y, peak1_pos,  peak2_x, peak2_y, peak2_pos, peak12_x, peak12_y)


        peak3_x, peak3_y, peak3_pos=select_near_max(spec_x,spec_y, 
                peak1_pos + 0.85 * peak1_3_mean_shift, peak1_pos + 1.15 * peak1_3_mean_shift,
                Calibration.FIT_RANGE_SMALL_MARGIN, Calibration.FIT_RANGE_BIG_MARGIN)
        high_e_result=PeakFitting.fit_asym_gaus(peak3_x,peak3_y,bin_width)
        #print(high_e_result.fit_report())
        peak_centers=np.array([low_e_1_result.params['cen'].value, low_e_2_result.params['cen'].value, high_e_result.params['cen'].value])
        #peak_center_err = np.array(
        #    [low_e_1_result.params['cen'].stderr, low_e_2_result.params['cen'].stderr, high_e_result.params['cen'].stderr])


        #peak_center_err[peak_center_err==None]=bin_width # set to 0.5 bin width if it is failed
        peak_center_err=np.array([0.25*bin_width, 0.25*bin_width, 0.25*bin_width])
        #print('line')
        #print(PHOTO_PEAKS_POS, peak_centers,peak_center_err)

        #peak_centers=np.array([
        #    np.sum(peak1_x*peak1_y)/np.sum(peak1_y),
        #    np.sum(peak2_x*peak2_y)/np.sum(peak2_y),
        #    np.sum(peak3_x*peak3_y)/np.sum(peak3_y)]) #weighted center, only select 

        #print(peak_centers)


        linfit_result=PeakFitting.fit_linear(PHOTO_PEAKS_POS, peak_centers,peak_center_err)
        #print(linfit_result.fit_report())
        try:
            slope_error=round(linfit_result.params['slope'].stderr,3)
            offset_error= round(linfit_result.params['intercept'].stderr, 4)
        except TypeError:
            print(f"ERROR: D {detector}-{pixel} slope, offset error not available.")
            slope_error=0.5*bin_width/(81-31)
            offset_error=0.5*bin_width
        best_fit={'slope': round(linfit_result.params['slope'].value,3),
                  'offset': round(linfit_result.params['intercept'].value,4),
                  'slope_error': slope_error,
                  'offset_error':offset_error
                  }
        if pdf:
            fig, axs=plt.subplots(1,2, figsize=(12,5))
            #axs[0].plot(peak1_x, low_e_1_result.best_fit, 'k',  label='31 keV peak fit')
            #axs[0].plot(peak2_x, low_e_2_result.best_fit, 'k',  label='35 keV peak fit')
            axs[0].plot(spec_x, spec_y, label='Obs')
            for _x in peak_centers:
                axs[0].axvline(_x, linewidth=1)
            #axs[0].plot(peak3_x, high_e_result.best_fit, 'k', label='81 keV peak fit')
            axs[0].set_xlabel('ADC channel')
            axs[0].set_ylabel('Counts')
            plt.suptitle(f'Detector {detector},pixel {pixel}, subspec {sbspec_id} ')
            axs[0].legend()
            axs[1].errorbar(PHOTO_PEAKS_POS, peak_centers, yerr=peak_center_err, fmt='o', label='Obs photopeaks')
            axs[1].plot(PHOTO_PEAKS_POS, linfit_result.best_fit, 'k',  label='Linear fit')
            text=f'baseline = {best_fit["offset"]:.3f} +/- {best_fit["offset_error"]}\n'\
                 f'gain = {best_fit["slope"]:.3f} +/- {best_fit["slope_error"]}\n'
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs[1].text(0.3, 0.95,
                     text,transform=axs[1].transAxes, fontsize=12,
                        verticalalignment='top',  bbox=props)

            axs[1].set_xlabel('Energy (keV)')
            axs[1].set_ylabel('Peak position')
            axs[1].legend()

            #plt.show()
            plt.tight_layout()

            pdf.savefig(fig)
            plt.close()
        return best_fit


    def process_one_run(self, calibration_id, create_pdf=True, pdf_path=DEFAULT_OUTPUT_DIR):
        docs = mdb.get_calibration_run_data(calibration_id)
        if not docs:
            print("Calibration run {} data doesn't exist".format(calibration_id))
            return
        data = docs[0]
        sbspec_formats = data['sbspec_formats']
        spectra = data['spectra']
        pdf=None
        if create_pdf:
            pdf_filename = os.path.abspath(
                os.path.join(Calibration.DEFAULT_OUTPUT_DIR, 'calibration_{}_gaus.pdf'.format(calibration_id)))
            pdf = PdfPages(pdf_filename)

        slope = np.zeros((32, 12))
        offset = np.zeros((32, 12))
        slope_error = np.zeros((32, 12))
        offset_error = np.zeros((32, 12))

        report = {}
        report['fit_parameters'] = []
        print('Processing calibration run {} ...'.format(calibration_id))

        for spec in spectra:
            # iterate over sub spectra
            detector, pixel, sbspec_id, start, bin_width, spectrum = spec
            print(f'{detector=},{pixel=},{sbspec_id=}')
            spectrum = np.array(spectrum)
            if np.sum(spectrum) < Calibration.PEAK_MIN_COUNTS:
                # total counts less than the limit
                continue

            end = start + bin_width * len(spectrum)
            if start > Calibration.FIT_RANGE[0] or end < Calibration.FIT_RANGE[1]:
                continue
            try:
                best_fit=self.find_solution(detector, pixel, sbspec_id, spectrum,start, bin_width,pdf)
            except:
                raise
                continue

            slope[detector][pixel] = best_fit['slope']
            offset[detector][pixel] = best_fit['offset']
            slope_error[detector][pixel] = best_fit['slope_error']
            offset_error[detector][pixel] = best_fit['offset_error']
        report['elut'] = self.compute_elut(offset, slope)
        slope1d = np.hstack(slope)
        offset1d = np.hstack(offset)
        slope_error_1d = np.hstack(slope_error)
        offset_error_1d = np.hstack(offset_error)

        # do calibration
        sum_spectra = {}#sum spectrum for all detectors
        # cc = TCanvas()
        xvals = []
        for spec in spectra:

            detector, pixel, sbspec_id, start, bin_width, spectrum = spec
            if sum(spectrum) < Calibration.PEAK_MIN_COUNTS:
                continue
            num_points = len(spectrum)

            spectrum=np.array(spectrum)

            end = start + bin_width * num_points  # end ADC

            if slope[detector][pixel] > 0 and offset[detector][pixel] > 0:
                #energies = (np.linspace(start, end - bin_width, num_points) - offset[detector][pixel]) / \
                #           slope[detector][pixel]
                spec_x, spec_y = self.get_subspec(spectrum, start, bin_width, False)
                energies=(spec_x-offset[detector][pixel])/slope[detector][pixel]
                if sbspec_id not in sum_spectra:
                    #first occurent of the spectrum, any detector
                    min_energy,max_energy=np.min(energies)*0.8,np.max(energies)*1.2
                    xvals = np.linspace(min_energy, max_energy,
                                        int((num_points + 1) * 1.4))
                    #xvals are the same for all detectors
                    sum_spectra[sbspec_id] = [[], []]  # sum spectrum x vs. y
                    sum_spectra[sbspec_id][0] = xvals
                    sum_spectra[sbspec_id][1] = np.zeros_like(xvals)

                yvals = interp(energies, spectrum / bin_width, xvals)
                sum_spectra[sbspec_id][1] += yvals

        sub_sum_spec = {}
        #produce sum spectrum
        points = 1150
        energy_range = np.linspace(-10, 450, points)
        sbspec_sum = np.zeros(points)
        for key, val in sum_spectra.items():  # mongodb doesn't support array
            sub_sum_spec['sbspec - {}'.format(key)] = [
                v.tolist() for v in sum_spectra[key]
            ]
            sbspec_sum += interp(sum_spectra[key][0], sum_spectra[key][1],
                                 energy_range)

        sub_sum_spec['sbspec sum'] = [energy_range.tolist(), sbspec_sum.tolist()]

        report['slope'] = slope1d.tolist()
        report['offset'] = offset1d.tolist()
        report['slope_error'] = slope_error_1d.tolist()
        report['offset_error'] = offset_error_1d.tolist()
        report['sum_spectra'] = sub_sum_spec
        if create_pdf:
            pdf.close()
            report['pdf']=pdf_filename
            print("PDF filename:",pdf_filename)
        # calibrated sum spectra


        mdb.update_calibration_analysis_report(calibration_id, report)


if __name__ == '__main__':
    #output_dir = DEFAULT_OUTPUT_DIR
    # output_dir='./'

    if len(sys.argv) == 1:
        print('Usage ./calibration_ecc <calibration_run_id> [<end id>]')
    elif len(sys.argv) >= 2:
        start_id = int(sys.argv[1])
        end_id = start_id
        p=Calibration()
        if len(sys.argv) >= 3:
            end_id = int(sys.argv[2])
        for i in range(start_id, end_id + 1):
            p.process_one_run(i)
