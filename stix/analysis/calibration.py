#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculation of energy calibration factors. 
This script relies on pyroot
It can be downloaded from http://root.cern.ch
As the pre-compiled version doesn't support python3, one needs to download the source code and compile on your local PC according to steps as below:
1. cmake 
cmake ../source   -Dpython3=ON -DPYTHON_EXECUTABLE=/usr/bin/python3 
-DPYTHON_INCLUDE_DIR=/usr/include/python3.8 -DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.8.so -DCMAKE_INSTALL_PREFIX=/opt/root6_python3/
2. make 
3. make install
"""

import os
import sys
import numpy as np
import time
from array import array
from datetime import datetime
from stix.core import stix_datatypes as sdt
from stix.core import mongo_db as db
from ROOT import TGraph, TFile, TCanvas, TH1F, gROOT, TBrowser, gSystem, TH2F, gPad, TF1, TGraphErrors, gStyle, TSpectrum, gRandom, TPaveLabel, TPaveText

from scipy.interpolate import interp1d

FIT_MIN_X = 252
FIT_MAX_X = 448
FIRST_PEAK_XMAX=350
MAX_ALLOWED_SIGMA_ERROR = 20  #maximum allowed peak error
MEAN_ENERGY_CONVERSION_FACTOR = 2.31
MIN_COUNTS = 100
#2.3 ADC/keV
#Estimated energy conversion factor

DEFAULT_OUTPUT_DIR = '/data/calibration/'
MIN_COUNTS_PEAK_FIND = 50
ELUT_ENERGIES = [
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36,
    40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150
]

#PHOTO_PEAKS_POS = [30.8, 35.2, 80.90] based spectral lines from olivier
PHOTO_PEAKS_POS = [30.85, 35.13, 81]

PRINT_TO_PDF = True

mdb = db.MongoDB()
gROOT.SetBatch(True)


def compute_elut(offset, slope):
    elut = []
    for det in range(0, 32):
        for pix in range(0, 12):
            p0 = offset[det][pix]
            p1 = slope[det][pix]
            #print(det, pix, p0, p1)

            if p0 > 0 and p1 > 0:
                row = [det, pix, p0, p1]
                Elows = [int(4 * (p0 + p1 * x)) for x in ELUT_ENERGIES]
                row.extend(Elows)
                elut.append(row)
    return elut





def interp(xvals, yvals, xnew):
    # x y define orignal points
    # xnew interpolated data points
    f2 = interp1d(xvals, yvals, bounds_error=False, fill_value=0)
    y = f2(xnew)
    return y


def create_graph_errors(x, y, ex, ey, title, xlabel="x", ylabel="y"):
    n = len(x)
    g = TGraphErrors(n, array('d', x), array('d', y), array('d', ex),
                     array('d', ey))
    g.GetXaxis().SetTitle(xlabel)
    g.GetYaxis().SetTitle(ylabel)
    g.SetTitle(title)
    return g


def create_graph(x, y, title="", xlabel="x", ylabel="y"):
    n = x.size
    g = TGraph(n, x.astype('float'), y.astype('float'))
    g.GetXaxis().SetTitle(xlabel)
    g.GetYaxis().SetTitle(ylabel)
    g.SetTitle(title)
    return g


def find_peaks(detector, pixel, subspec, start, num_summed, spectrum, fo):
    gStyle.SetOptFit(111)
    num_points=spectrum.size
    x_full=np.linspace(start, start+num_points*num_summed, 
            num_points)+0.5*num_summed
    sel=(x_full> FIT_MIN_X) & (x_full<FIT_MAX_X) 

    x,y = x_full[sel], spectrum[sel]

    #spectrum in the predefined range
    if x.size==0:
        print('Can not find sub spectrum of ERROR:', detector, pixel)
        return None, None

    total_counts = np.sum(y)
    if total_counts < MIN_COUNTS:
        print('Too less counts:', detector, pixel)
        return None, None

    name = '{}_{}_{}'.format(detector, pixel, subspec)
    title = 'detector {} pixel {} subspec {}'.format(detector, pixel, subspec)
    g_full_spec = create_graph(x_full, spectrum,
                         'Original spec - {}'.format(name), 'ADC channel',
                         'Counts')
    peak1_y = np.max(y[x<FIRST_PEAK_XMAX])
    peak1_x = x[y==peak1_y][0]
    #find the peak with highest counts in the range

    x_shift = MEAN_ENERGY_CONVERSION_FACTOR * (81 - 31)

    # max conversion factor = 2.5 ADC/keV
    fit_range_x_left = 5
    fit_range_x_right = 15
    fit_range_peak3_x_left = 2


    peak3_sel = (x > peak1_x + 0.9 * x_shift) & (x<peak1_x + 1.1 * x_shift)
    peak3_x=x[peak3_sel] 
    peak3_y=y[peak3_sel] 
    peak3_max_y = np.max(peak3_y)
    peak3_max_x = peak3_x[peak3_y==peak3_max_y][0]

    peak2_x = peak1_x + 4.2 * MEAN_ENERGY_CONVERSION_FACTOR


    fgaus1 = TF1('fgaus1_{}'.format(name), 'gaus', peak1_x - fit_range_x_left,
                 peak1_x + fit_range_x_right)
    fgaus2 = TF1('fgaus2_{}'.format(name), 'gaus', peak2_x - fit_range_x_left,
                 peak2_x + fit_range_x_right)
    fgaus3 = TF1('fgaus3_{}'.format(name), 'gaus',
                 peak3_max_x - fit_range_peak3_x_left,
                 peak3_max_x + fit_range_x_right)

    gspec = create_graph(x, y, 'Spectrum - {}'.format(title), 'ADC channel',
                   'Counts')

    gspec.Fit(fgaus1, 'RQ')
    gspec.Fit(fgaus2, 'RQ+')
    gspec.Fit(fgaus3, 'RQ')
    par1 = fgaus1.GetParameters()
    par2 = fgaus2.GetParameters()
    par3 = fgaus3.GetParameters()
    par3_errors = fgaus3.GetParErrors()

    fgaus12 = TF1('fgaus12_{}'.format(name), 'gaus(0)+gaus(3)', par1[1] - 2,
                  par2[1] + 3, 6)
    par = array('d', [par1[0], par1[1], par1[2], par2[0], par2[1], par2[2]])
    fgaus12.SetParameters(par)
    gspec.Fit(fgaus12, 'RQ+')
    param = fgaus12.GetParameters()
    param_errors = fgaus12.GetParErrors()

    fo.cd()
    gspec.Write(f'spec_fits_{name}')
    g_full_spec.Write(f'spec_{name}')
    fgaus12.Write()
    result = {'detector': detector, 'pixel': pixel, 'sbspec_id': subspec}
    try:
        result['peaks'] = {
            'peak1': (param[0], param[1], param[2]),
            'peak2': (param[3], param[4], param[5]),
            'peak3': (par3[0], par3[1], par3[2]),
            'peak1error': (param_errors[0], param_errors[1], param_errors[2]),
            'peak2error': (param_errors[3], param_errors[4], param_errors[5]),
            'peak3error': (par3_errors[0], par3_errors[1], par3_errors[2])
        }
    except Exception as e:
        print(str(e))
    peak_x = []
    peak_ex = []
    peak_y = []
    peak_ey = []
    if param_errors[2] < MAX_ALLOWED_SIGMA_ERROR:
        peak_x.append(PHOTO_PEAKS_POS[0])
        peak_ex.append(0.)
        peak_y.append(param[1])
        peak_ey.append(param_errors[1])
    if par3_errors[2] < MAX_ALLOWED_SIGMA_ERROR:
        #compensation=0.5 #we observed that there is about 0.5 adc channels shift if it  is fitting with single gaussian function
        compensation=0
        peak_x.append(PHOTO_PEAKS_POS[2])
        peak_ex.append(0.125)
        peak_y.append(par3[1]+compensation)
        peak_ey.append(par3_errors[1])
    #peak_x=[30.8, 34.9, 81]
    #peak_ex=[.0, 0., 0.]
    gpeaks = None
    if len(peak_x) >= 2:
        gpeaks = create_graph_errors(peak_x, peak_y, peak_ex, peak_ey, title,
                              'Energy (keV)', 'Peak position (ADC)')
        gpeaks.Fit('pol1', 'Q')
        gpeaks.GetYaxis().SetRangeUser(0.9 * peak_y[0], peak_y[-1] * 1.1)
        gpeaks.Write('gpeaks_{}'.format(name))

        calibration_params = gpeaks.GetFunction('pol1').GetParameters()
        chisquare = gpeaks.GetFunction('pol1').GetChisquare()
        fcal_errors = gpeaks.GetFunction('pol1').GetParErrors()
        result['fcal'] = {
            'p0': calibration_params[0],
            'p1': calibration_params[1],
            'chi2': chisquare,
            'p0error': round(fcal_errors[0], 5),
            'p1error': round(fcal_errors[1], 5),
        }

    return result, [g_full_spec, gspec, gpeaks, fgaus12, fgaus3]


def process_one_run(calibration_id, create_pdf=True, pdf_path=DEFAULT_OUTPUT_DIR):
    runs = list(mdb.get_calibration_run_data(calibration_id))
    if not runs:
        print("Calibration run {} doesn't exist".format(calibration_id))
        return
    data=runs[0]

    sbspec_formats = data['sbspec_formats']
    spectra = data['spectra']

    fname_out = os.path.abspath(
        os.path.join(pdf_path, 'calibration_{}'.format(calibration_id)))

    f = TFile("{}.root".format(fname_out), "recreate")

    is_top = True

    slope = np.zeros((32, 12))
    offset = np.zeros((32, 12))
    slope_error = np.zeros((32, 12))
    offset_error = np.zeros((32, 12))

    report = {}
    report['fit_parameters'] = []
    print('Processing calibration run {} ...'.format(calibration_id))

    canvas = TCanvas("c", "canvas", 1200, 500)
    pdf = '{}.pdf'.format(fname_out)
    if PRINT_TO_PDF:
        print(f'Plots will be written to {pdf}')
        canvas.Print(pdf + '[')
    # make cover
    cover = TCanvas()
    t1 = TPaveText(.1, .6, .9, .9)
    t1.AddText(f"Calibration run {calibration_id} analysis report")
    t1.SetTextAlign(22)
    t1.SetTextFont(52)
    t1.SetTextColor(4)
    t1.SetFillColor(24)
    t1.Draw()
    t2ptxt = TPaveText(.1, .3, .9, .58)
    t2ptxt.SetTextAlign(12)
    t2ptxt.SetTextFont(52)
    now = f'Created at {datetime.now().isoformat()}'
    t2ptxt.AddText(now)
    t2ptxt.Draw()
    if PRINT_TO_PDF:
        cover.Print(pdf)

    canvas.Divide(3, 2)
    last_plots = None
    for spec in spectra:
        #spectrum data from database
        if sum(spec[5]) < MIN_COUNTS_PEAK_FIND:
            continue
        detector = spec[0]
        pixel = spec[1]
        sbspec_id = spec[2]
        start = spec[3]
        num_summed = spec[4]
        end = start + num_summed * len(spec[5])
        spectrum = np.array(spec[5])
        if start > FIT_MAX_X or end < FIT_MIN_X:
            #break
            continue
        par, plots = find_peaks(detector, pixel, sbspec_id, start, num_summed,
                                spectrum, f)
        if not par and not plots:
            continue
        if last_plots:
            canvas.cd(1)
            if last_plots[0]:
                last_plots[0].Draw("AL")
            canvas.cd(2)
            if last_plots[1]:
                last_plots[1].Draw("AL")
            canvas.cd(3)
            if last_plots[2]:
                last_plots[2].Draw("AP")
            canvas.cd(4)
            if plots[0]:
                plots[0].Draw("AL")
            canvas.cd(5)
            if plots[1]:
                #plot[3].Draw()
                plots[1].Draw("AL")
            canvas.cd(6)
            if plots[2]:
                plots[2].Draw("AP")
            if PRINT_TO_PDF:
                canvas.Print(pdf)
            last_plots = []
        else:
            last_plots = plots
        report['fit_parameters'].append(par)

        if par:
            if 'fcal' in par:
                slope[detector][pixel] = par['fcal']['p1']
                offset[detector][pixel] = par['fcal']['p0']
                slope_error[detector][pixel] = par['fcal']['p1error']
                offset_error[detector][pixel] = par['fcal']['p0error']

    report['pdf'] = pdf
    report['elut'] = compute_elut(offset, slope)

    slope1d = slope.flatten()
    offset1d = offset.flatten()
    slope_error_1d = slope_error.flatten()
    offset_error_1d = offset_error.flatten()

    #do calibration
    sum_spectra = {}
    cc = TCanvas()
    xvals = []
    for spec in spectra:
        if sum(spec[5]) < MIN_COUNTS_PEAK_FIND:
            continue
        detector = spec[0]
        pixel = spec[1]
        sbspec_id = spec[2]
        spectrum = np.array(spec[5])
        num_points = spectrum.size

        start = spec[3] #ADC channel start

        num_summed = spec[4]

        end = start + num_summed * num_points #end ADC 

        print(detector, pixel)
        if slope[detector][pixel] > 0 and offset[detector][pixel] > 0:
            energies = (np.linspace(start, end - num_summed, num_points) - offset[detector][pixel]) / slope[detector][pixel]
            if sbspec_id not in sum_spectra:
                min_energy = (start - offset[detector][pixel]
                              ) / slope[detector][pixel] * 0.8 #20% margin
                max_energy = (end - offset[detector][pixel]
                              ) / slope[detector][pixel] * 1.2 #20% margin
                xvals = np.linspace(min_energy, max_energy,
                                    int((num_points + 1) * 1.4)) 
                sum_spectra[sbspec_id] = [[], []] #sum spectrum x vs. y
                sum_spectra[sbspec_id][0] = xvals
                sum_spectra[sbspec_id][1] = np.zeros(len(xvals))

            yvals = interp(energies, spectrum / num_summed, xvals)
            sum_spectra[sbspec_id][1] += yvals

       
    sub_sum_spec = {}

    points = 1150
    energy_range = np.linspace(-10, 450, points)
    sbspec_sum = np.zeros(points)
    for key, val in sum_spectra.items():  #mongodb doesn't support array
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
    #calibrated sum spectra

    mdb.update_calibration_analysis_report(calibration_id, report)

    hist_slope = TH1F(
        "hist_slope",
        "Energy conversion factors; Conversion factors (ADC / keV); Counts",
        100, 0.8 * np.min(slope1d), 1.2 * np.max(slope1d))
    for s in slope1d:
        hist_slope.Fill(s)
    hist_offset = TH1F("hist_offset", "Baseline; Baseline (ADC); Counts", 100,
                       0.8 * np.min(offset1d), 1.2 * np.max(offset1d))
    for s in offset1d:
        hist_offset.Fill(s)
    ids = np.arange(384)
    g_slope = create_graph(ids, slope1d, 'conversion factor', ' pixel #',
                     'conversion factor')
    g_offset = create_graph(ids, offset1d, 'baseline', ' pixel #', 'baseline')
    c2 = TCanvas()
    c2.Divide(2, 2)
    c2.cd(1)
    hist_slope.Draw('hist')
    c2.cd(2)
    hist_offset.Draw('hist')
    c2.cd(3)
    g_slope.Draw('AL')
    c2.cd(4)
    g_offset.Draw('AL')
    if PRINT_TO_PDF:
        c2.Print(pdf)
        canvas.Print(pdf + ']')

    hist_slope.Write("hist_slope")
    hist_offset.Write("hist_offset")
    g_slope.Write("g_slope")
    g_offset.Write("g_offset")
    if PRINT_TO_PDF:
        print('done.\nFile {} generated'.format(pdf))

    f.Close()




if __name__ == '__main__':
    pdf_path= DEFAULT_OUTPUT_DIR
    #output_dir='./'
    if len(sys.argv) == 1:
        print("Usage ./calibration <run_id>")
    elif len(sys.argv) >= 2:
        start_id = int(sys.argv[1])
        end_id=start_id
        if len(sys.argv) >= 3:
            end_id = int(sys.argv[2])
        for i in range(start_id, end_id + 1):
            process_one_run(i, pdf_path)
