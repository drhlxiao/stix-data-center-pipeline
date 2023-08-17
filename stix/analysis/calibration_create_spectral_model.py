#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

from matplotlib.backends.backend_pdf import PdfPages
import os
import sys

import numba
import numpy as np
import time
from array import array

from scipy.interpolate import interp1d
from pathlib import Path
from matplotlib import pyplot as plt
from scipy import stats

from stix.core import datatypes as sdt
from stix.core import mongo_db as db
from stix.analysis import config
import pickle

mdb = db.MongoDB()

#
ELUT_TIME = '2021-08-29T00:00:00'

# 2.3 ADC/keV
# Estimated energy conversion factor
# import numba
SPECTRUM_MODEL_DATA_FILES = Path(
    __file__).parent.parent / 'data' / 'ExpSpecModel.npz'


def interp(xvals, yvals, xnew):
    # x y define orignal points
    # xnew interpolated data points
    f2 = interp1d(xvals, yvals, bounds_error=False, fill_value=0)
    y = f2(xnew)
    return y


def chi2test(obs_y: np.ndarray, exp_y: np.ndarray):
    """
        Do chisquare test between to spectra
    Arguments
    obs_y: numpy.ndarray
        observed energy spectrum
    exp_y: numpy.ndarray
        expected energy spectrum
    Returns
        chisquare: float
            chiquare
        dof:
            dof: dregress of freedom
    """
    # exp_y: expected spectrum
    # obs_y: observed spectra
    norm_factor = np.sum(obs_y) / np.sum(exp_y)
    exp_y[exp_y == 0] = 1e-12
    dof = obs_y.size - 1
    return np.sum((obs_y - exp_y * norm_factor)**2 / exp_y), dof
    # stats.chisquare is too slow
    # return stats.chisquare(f_obs=obs_y, f_exp=exp_y)


class Calibration(object):
    """
    Fit calibration spectra using  spectrum models loaded from npy file

    """
    FIT_RANGE_KEV = (28, 85)  # 20 keV - 90 keV
    SPECTRUM_SEL_ADC_RANGE = (252, 448)
    # used to check if do fitting for a subspectrum
    OFFSETS_LIMITS = (210, 235)  # default range, load from ELUT by default
    MIN_COUNTS = 100
    DEFAULT_OUTPUT_DIR = '/data/calibration/'
    STEPS = 2000  # number of data points generated for each iteration
    MAX_DEPTH = 1  # random grid recursive search depth
    SCALE_FACTOR = 10  # shrink search region after each iteration
    PEAK_MIN_COUNTS = 50  # spectrum with little counts ignored
    MEAN_GAIN = 2.31
    ELUT_ENERGIES = [
        4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32,
        36, 40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150
    ]
    SUM_SPECTRUM_STEPS = 10000
    PIXEL_EBINS = np.linspace(0, 100, SUM_SPECTRUM_STEPS)

    def __init__(self):
        # self.pixel_fspec_models, self.default_slopes, self.default_offsets = \
        #    self.load_spectrum_models(SPECTRUM_MODEL_DATA_FILES)
        pass

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
        max_count = np.max(spectrum)
        peak_x = _x[spectrum == max_count][0]
        min_adc = peak_x - (
            31 - Calibration.FIT_RANGE_KEV[0]) * Calibration.MEAN_GAIN
        max_adc = peak_x + (
            -31 + Calibration.FIT_RANGE_KEV[1]) * Calibration.MEAN_GAIN

        if only_fit_range:
            clip = (_x > min_adc) & (_x < max_adc)
            return _x[clip], spectrum[clip]
        return _x, spectrum

    def get_intensity_from_model(self, detector: int, pixel: int,
                                 offset: float, slope: float,
                                 spec_adcs: np.ndarray) -> object:
        """
        convert adc channels to keV
        Args
            offset: float
                offset
            slope: float
                gain
        returns
            an array having the same length as adc spectrum
        """
        energies = (spec_adcs - offset) / slope
        return energies, self.pixel_fspec_models[detector * 12 +
                                                 pixel](energies)

    def create_sum_spectra(self, calibration_id, spectra, slope: np.ndarray,
                           offset: np.ndarray):
        # do calibration
        sum_spectra = {}  # sum spectrum for all detectors
        # cc = TCanvas()

        pixel_sum_spectra = np.zeros((32, 12, Calibration.SUM_SPECTRUM_STEPS
                                      ))  # pixel sum spectrum, 100 energy bins

        xvals = []
        for spec in spectra:
            print(spec)

            detector, pixel, sbspec_id, start, bin_width, spectrum = spec
            num_points = len(spectrum)
            spectrum = np.array(spectrum)
            if np.sum(spectrum) < Calibration.PEAK_MIN_COUNTS:
                continue
            end = start + bin_width * num_points  # end ADC
            if slope[detector][pixel] > 0 and offset[detector][pixel] > 0:
                # energies = (np.linspace(start, end - bin_width, num_points) - offset[detector][pixel]) / \
                #           slope[detector][pixel]
                spec_x, spec_y = self.get_subspec(spectrum, start, bin_width,
                                                  False)
                energies = (spec_x -
                            offset[detector][pixel]) / slope[detector][pixel]
                if sbspec_id not in sum_spectra:
                    # first occurent of the spectrum, any detector
                    min_energy, max_energy = np.min(energies) * 0.8, np.max(
                        energies) * 1.2
                    xvals = np.linspace(min_energy, max_energy,
                                        int((num_points + 1) * 1.4))
                    sum_spectra[sbspec_id] = [xvals,
                                              np.zeros_like(xvals)
                                              ]  # sum spectrum x vs. y
                yvals = interp(energies, spectrum / bin_width, xvals)
                sum_spectra[sbspec_id][1] += yvals

                amps = interp(energies, spectrum / bin_width,
                              Calibration.PIXEL_EBINS)
                pixel_sum_spectra[detector][pixel] += amps
        sub_sum_spec = {}
        points = 1150
        energy_range = np.linspace(-10, 450, points)
        sbspec_sum = np.zeros(points)
        for key, val in sum_spectra.items():  # mongodb doesn't support array
            sub_sum_spec[f'sbspec - {key}'] = [
                v.tolist() for v in sum_spectra[key]
            ]
            sbspec_sum += interp(sum_spectra[key][0], sum_spectra[key][1],
                                 energy_range)

        sub_sum_spec['sbspec sum'] = [
            energy_range.tolist(), sbspec_sum.tolist()
        ]
        return sub_sum_spec, pixel_sum_spectra

    def process_one_run(self, calibration_id, create_pdf=True):
        docs = mdb.get_calibration_run_data(calibration_id)
        if not docs:
            print(
                "Calibration run {} data doesn't exist".format(calibration_id))
            return
        data = docs[0]
        spectra = data['spectra']
        pdf, pdf_filename = None, None

        if create_pdf:
            pdf_filename = os.path.abspath(
                os.path.join(Calibration.DEFAULT_OUTPUT_DIR,
                             f'calibration_{calibration_id}_chisquare.pdf'))
            pdf = PdfPages(pdf_filename)
        slope = np.zeros((32, 12))
        offset = np.zeros((32, 12))
        chi2 = np.zeros((32, 12))
        slope_error = np.zeros((32, 12))
        offset_error = np.zeros((32, 12))
        report = {}
        report['fit_parameters'] = []
        print(f'Processing calibration run {calibration_id} ...')

        # from ELUT
        # elut=config.Elut.get_onboard_elut(ELUT_TIME)
        # slope=np.array(elut['slopes']).reshape((32,12))
        # offset=np.array(elut['offsets']).reshape((32,12))

        slope = np.array(data['analysis_report']['slope']).reshape((32, 12))
        offset = np.array(data['analysis_report']['offset']).reshape((32, 12))
        sub_sum_spec, pixel_sum_spectra = self.create_sum_spectra(
            calibration_id, spectra, slope, offset)

        calibration_result_filename = os.path.join(
            Calibration.DEFAULT_OUTPUT_DIR,
            f'calibration_results_{calibration_id}.npz')
        print(f'calibration_results saved to: {calibration_result_filename}')

        np.savez(calibration_result_filename, pixel_sum_spectra, slope, offset)


if __name__ == '__main__':
    # output_dir = DEFAULT_OUTPUT_DIR
    # output_dir='./'
    if len(sys.argv) == 1:
        print('Usage ./calibration_ecc <calibration_run_id> [<end id>]')
    elif len(sys.argv) >= 2:
        cal = Calibration()
        start_id = int(sys.argv[1])
        end_id = start_id
        if len(sys.argv) >= 3:
            end_id = int(sys.argv[2])
        for i in range(start_id, end_id + 1):
            cal.process_one_run(i)
