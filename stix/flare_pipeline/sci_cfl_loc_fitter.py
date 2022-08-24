#!/usr/bin/python3
'''

This script should run after running flare_sci_L1_analysis because it requires data prepared by flare sci


L1 flare data analysis pipeline
  what it does here includes
  - do background subtraction for flare data
  - find flare location solution
  - pre-requirements:
     flare data should be processed by flare_goes_class.py and flare_detection.py and sci_packets_analyzer.py 
    1) flare_goes_class prepare ephemeris for flare data
    2) sci_packets_analyzer prepares spectra and  check if it is background
    3) flare detection tells flare start time, end time, flare peak counts
Author: Hualin Xiao, Aug. 2021
email: hualin.xiao@fhnw.ch
'''
import sys
import os
import json
import numpy as np
from datetime import datetime
from stix.core import datatypes as sdt
from stix.spice import datetime as st
from stix.core import mongo_db as db
from stix.core import logger
from stix.core import config
from stix.flare_pipeline import flare_location_fitter as flf
from pprint import pprint
from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Rectangle, PathPatch
import matplotlib.colors as colors
from matplotlib.path import Path
from stix.utils import bson
from stix.utils import energy_bins as seb
from stix.flare_pipeline import flare_spice as fsp
import matplotlib

matplotlib.use('Agg')

DET_ID_START = 20
DET_ID_END = 32  #only detector 20-32 selected  for detector mean counts calculation
DET_SENSITVE_AREA = 80.96  #detector pixel sensitive area in units of mm2
BKG_MIN_DURATION = 120  #minimal duration of bsd that can be used  for bkg subtraction, to reduce systematic error
GRID_OPEN_AREA_RATIO = np.array([
    0.383332095, 0.243715836, 0.260355706, 0.263845616, 0.284908432,
    0.248516311, 0.411123486, 0.256797947, 1, 1, 0.697247857, 0.308301538,
    0.335181618, 0.263335304, 0.244122951, 0.261332516, 0.308298271,
    0.335179274, 0.308303277, 0.266539578, 0.254592587, 0.256522323,
    0.280037904, 0.239474758, 0.291241432, 0.256524112, 0.231796481,
    0.271006556, 0.390570976, 0.260070316, 0.264345285, 0.261086447
])
#// = slit_width_front/pitch_front*slit_width_rear/pitch_rear
#		// nominal parameters are used for gird 12,16,17,18 , real measurements are used for the rest

flare_pipeline_out_path = config.get_config(
    'pipeline.daemon.flare_pipeline_path')
mdb = db.MongoDB()


class FlareDataAnalyzer(object):

    def __init__(self):
        self.bsd_db = mdb.get_collection_bsd()
        self.flare_db = mdb.get_collection("flares")

    def get_background_request_data(self,
                                    signal_unix_time,
                                    emin,
                                    emax,
                                    background_bsd_id=None):
        #load a bulk science data report for background subtraction or find the closest in time L1 data for background subtraction
        if background_bsd_id is not None:
            return self.bsd_db.find_one({'_id': background_bsd_id})
        '''
        find background data by database query
        bsd_ealier=list(self.bsd_db.find({'synopsis.is_background':True,
            'synopsis.duration':{'$gt':BKG_MIN_DURATION},
            'synopsis.emin':{'$lte':emin},
            'synopsis.emax':{'$gte':emax},
            'start_unix':{'$lt':signal_unix_time}}).sort('start_unix',-1).limit(1))

        bsd_later=list(self.bsd_db.find({'synopsis.is_background':True, 
            'synopsis.duration':{'$gt':BKG_MIN_DURATION},
            'synopsis.emin':{'$lte':emin},
            'synopsis.emax':{'$gte':emax},
            'start_unix':{'$gt':signal_unix_time}}).sort('start_unix',1).limit(1))
      '''
        #find background data by database query
        bsd_ealier = list(
            self.bsd_db.find({
                'request_form.purpose':
                'Background',  #only select those marked as background
                'synopsis.duration': {
                    '$gt': BKG_MIN_DURATION
                },
                'synopsis.emin': {
                    '$lte': emin
                },
                'name': 'L1',
                'synopsis.emax': {
                    '$gte': emax
                },
                'start_unix': {
                    '$lt': signal_unix_time
                }
            }).sort('start_unix', -1).limit(1))

        bsd_later = list(
            self.bsd_db.find({
                'request_form.purpose':
                'Background',  #only select those marked as background
                'synopsis.duration': {
                    '$gt': BKG_MIN_DURATION
                },
                'name': 'L1',
                'synopsis.emin': {
                    '$lte': emin
                },
                'synopsis.emax': {
                    '$gte': emax
                },
                'start_unix': {
                    '$gt': signal_unix_time
                }
            }).sort('start_unix', 1).limit(1))
        try:
            t0 = np.abs(bsd_ealier[0]['start_unix'] - signal_unix_time)

        except (KeyError, IndexError):
            t0 = np.inf
        try:
            t1 = np.abs(bsd_later[0]['start_unix'] - signal_unix_time)
        except (KeyError, IndexError):
            t1 = np.inf
        closet = min(t1, t0)

        if closet == np.inf:
            return None
        elif closet == t0:
            return bsd_ealier[0]
        else:
            return bsd_later[0]

    def process_one_flare(self,
                          bsd_doc,
                          bkg_bsd_id=None,
                          plot_fig=True,
                          save_fig=True):
        """tasks:
        1) background subtraction using pre-processed data 
        2) calculate CFL location 
        3) earth viewing angle of the flare
           """

        if not bsd_doc:
            print('Not a valid STIX L1 sci dataset')
        results = {}
        _id = bsd_doc['_id']
        print(f"Processing flare {_id}")
        try:
            flare_ids = bsd_doc['synopsis']['flares']['flare_ids']
        except KeyError:
            print(
                f'no flare information found in BSD #{_id} or it has not been processed by sci_l1_analyzer.py'
            )
            return None

        if len(flare_ids) > 1:
            print(
                f'BSD {_id} contains several flares,but only the first will be processed!'
            )
        flare_id = flare_ids[0]
        signal_energy_integrated_counts = np.array(
            bsd_doc['synopsis']['flares']['flare_pixel_counts']['counts'])
        #already integrated over energies,  1d vector 384 elements

        flare_emin = bsd_doc['synopsis']['flares']['flare_pixel_counts'][
            'emin']  #energy range selected for flare location calculation
        flare_emax = bsd_doc['synopsis']['flares']['flare_pixel_counts'][
            'emax']
        signal_duration = bsd_doc['synopsis']['flares']['integration_times'][0]
        bkg_doc = self.get_background_request_data(bsd_doc['start_unix'],
                                                   flare_emin, flare_emax,
                                                   bkg_bsd_id)

        if not bkg_doc:
            print(f'No background found for {_id}')
            return
        print(f'Using background: bsd #{bkg_doc["_id"]}')

        sig_emin, sig_emax = bsd_doc['synopsis']['emin'], bsd_doc['synopsis'][
            'emax']
        bkg_emin, bkg_emax = bkg_doc['synopsis']['emin'], bkg_doc['synopsis'][
            'emax']

        bkg_start_datetime = st.unix2datetime(bkg_doc['start_unix'])
        bkg_start_date = bkg_start_datetime.strftime('%Y-%m-%d')

        bkg_rates_full_energy = np.array(
            bkg_doc['synopsis']['pixel_rate_spec']
        )  #rate are not energy integrated, dimemsions 384x32
        bkg_duration = bkg_doc['synopsis']['duration']

        #calculate duration
        signal_rates_full_energy = np.array(
            bsd_doc['synopsis']['pixel_rate_spec']
        )  #rate are not energy integrated, dimemsions 384x32
        bsd_duration = bkg_doc['synopsis']['duration']

        e_range = [max(sig_emin, bkg_emin), min(sig_emax, bkg_emax)]

        bkg_spectrum = np.sum(bkg_rates_full_energy, axis=0) * bsd_duration

        sig_spectrum = np.sum(signal_rates_full_energy, axis=0) * bsd_duration

        spectrogram_data = bsd_doc['synopsis'][
            'spectrogram']  #rate are not energy integrated, dimemsions 384x32

        bkg_subtracted_spectrum = np.sum(
            signal_rates_full_energy - bkg_rates_full_energy,
            axis=0) * bsd_duration
        bkg_subtracted_spectrum_error = np.sqrt(
            np.sum(signal_rates_full_energy * bsd_duration +
                   bkg_rates_full_energy * bsd_duration,
                   axis=0))

        bkg_subtracted_spectrum[:e_range[0]] = 0
        bkg_subtracted_spectrum[e_range[1]:] = 0
        #the subtraction is not valid for e>emax and e<emin
        #print(f'{bkg_subtracted_spectrum=}')

        bkg_energy_integrated_rates = np.sum(
            bkg_rates_full_energy[:, flare_emin:flare_emax], axis=1
        )  #only select bins in the energy range, energy integrated counts

        print("bkg shape:", bkg_rates_full_energy.shape)
        bkg_rate_error = np.sqrt(
            bkg_energy_integrated_rates * bkg_duration) / bkg_duration
        # rate*duration =  original counts;   only statistic error considered

        bkg_cfl_rate = bkg_energy_integrated_rates[8 * 12:9 * 12]
        #print(f'{bkg_cfl_rate=}')

        bkg_cfl_rate_error = bkg_rate_error[8 * 12:9 * 12]
        #print(f'{bkg_cfl_rate_error=}')

        cfl_pixel_counts = signal_energy_integrated_counts[
            8 * 12:9 * 12]  #already energy bin integrated
        #print(f'{cfl_pixel_counts=}')

        cfl_pixel_bkgsub_counts = cfl_pixel_counts - signal_duration * bkg_cfl_rate
        #sqrt(cfl_counts)**2 +signal_duration* bkg_cfl_rate_error **2

        cfl_pixel_bkgsub_counts[cfl_pixel_bkgsub_counts < 0] = 0
        #reset negative counts to zeros, negative counts do not make physically

        cfl_pixel_bkgsub_count_error = np.sqrt(cfl_pixel_counts +
                                               signal_duration *
                                               bkg_cfl_rate_error**2)

        sig_bkgsub_counts = signal_energy_integrated_counts - signal_duration * bkg_energy_integrated_rates  #background subtracted counts

        sig_bkgsub_counts_error = np.sqrt(signal_energy_integrated_counts +
                                          signal_duration * bkg_rate_error**2)

        #print(f'sig shape, {sig_bkgsub_counts.shape}')

        counts_to_fluence_factors = 1 / (GRID_OPEN_AREA_RATIO *
                                         DET_SENSITVE_AREA)

        det_bkgsub_counts = np.array([
            np.sum(sig_bkgsub_counts[idet * 12:(idet * 12 + 12)])
            for idet in range(0, 32)
        ])
        det_bkgsub_counts_error = np.array([
            np.sqrt(
                np.sum(sig_bkgsub_counts_error[idet * 12:(idet * 12 + 12)]**2))
            for idet in range(0, 32)
        ])

        det_bkgsub_fluence, det_bkgsub_fluence_error = det_bkgsub_counts * counts_to_fluence_factors, det_bkgsub_counts_error * counts_to_fluence_factors
        #detector counts errors, only consider statistical uncertainties

        mean_fluence = np.mean(det_bkgsub_fluence[DET_ID_START:DET_ID_END])
        mean_fluence_error = np.sqrt(
            np.sum([
                det_bkgsub_fluence_error[i]**2
                for i in range(DET_ID_START, DET_ID_END)
            ]))
        flare_utc = st.unix2utc(bsd_doc['start_unix'])
        #print('start fitting')

        chi2_map, res = flf.fit_location(cfl_pixel_bkgsub_counts,
                                         cfl_pixel_bkgsub_count_error,
                                         mean_fluence,
                                         mean_fluence_error,
                                         flare_utc,
                                         use_small_pixels=True,
                                         use_det_fluence=True)

        xloc = res['min_pos'][0]
        yloc = res['min_pos'][1]
        print(f"Flare location:({xloc}, {yloc})")
        flare_spice_data = fsp.get_flare_spice(xloc,
                                               yloc,
                                               flare_utc,
                                               observer='solo')

        emission_earth = flare_spice_data['theta_flare_norm_earth_deg']
        emission_solo = flare_spice_data['theta_flare_norm_solo_deg']
        fig_filename = os.path.join(
            flare_pipeline_out_path,
            f'bsd_{_id}_flare_{flare_id}_l1_quick_analysis.png')

        result = {
            #'solution': res, # raise document too large exception
            'background_bsd_id': bkg_doc['_id'],
            'bkg_start_unix': bkg_doc['start_unix'],
            'fig_filename': fig_filename,
            'cfl_xloc': xloc,
            'cfl_yloc': yloc,
            'earth_flare_norm_deg': emission_earth,
            'earth_flare_solo_deg': flare_spice_data['earth_flare_solo_deg'],
            'solo_flare_norm_deg': emission_solo,
            'spice': flare_spice_data
        }
        result_doc = bson.dict_to_json(result)

        #store  results in database
        #print("Updating database...")
        self.bsd_db.update_one({'_id': bsd_doc['_id']},
                               {'$set': {
                                   'synopsis.flares.cfl': result_doc
                               }})

        #rate are not energy integrated, dimemsions 384x32
        result_doc['signal_bsd_id'] = _id
        self.flare_db.update_many({'flare_id': flare_id},
                                  {'$set': {
                                      'L1_pipeline': result_doc
                                  }})
        #print(f'{flare_id=}')
        """
        if plot_fig or save_fig:
            ds='steps-mid'
            fig, axs = plt.subplots(2,2, figsize=(13,10))
            
            spec_timebins=np.array([datetime.unix2datetime(xx) for xx in spectrogram_data['timebins']])
            spec_ebins=np.array([ seb.to_keV(x[0], x[1]) for x in spectrogram_data['ebins']])
            spec_data=np.transpose(np.array(spectrogram_data['data']))
            flare_start=datetime.unix2datetime(bsd_doc['synopsis']['flares']['start_unix_times'][0])
            flare_end=datetime.unix2datetime(bsd_doc['synopsis']['flares']['end_unix_times'][0])
            XX, YY=np.meshgrid(spec_timebins, spec_ebins)
            spec_data[spec_data<=0]=0.00001 #set to zero, negative dosen't make sense but could appear when doing bkg subtraction
            norm=colors.LogNorm(vmin=np.min(spec_data), vmax=np.max(spec_data))
            try:
                im2=axs[0,0].pcolormesh(XX,YY, spec_data, cmap='plasma', norm=norm)
            except Exception as e:
                print(e)
                im2=axs[0,0].imshow(spec_data)

            #fig.colorbar(im2, ax=axs[0,0])
            
            axs[0,0].set_title(f'Raw spectrogram (BSD #{bsd_doc["_id"]})')
            if flare_start < spec_timebins[0]:
                flare_start=spec_timebins[0]
            if flare_end> spec_timebins[-1]:
                flare_end=spec_timebins[-1]

            axs[0,0].vlines([flare_start, flare_end], ymin=flare_emin, ymax=flare_emax, colors='cyan')
            axs[0,0].set_yticklabels(spec_ebins)
            axs[0,0].set_xlabel(f'Start time {spec_timebins[0].strftime("%Y-%m-%dT%H:%M%S")}')



            axs[0,1].errorbar(range(32),
                    bkg_subtracted_spectrum, yerr=bkg_subtracted_spectrum_error,ds=ds, label='after bkg. sub.')
            axs[0,1].errorbar(range(32), bkg_spectrum  , yerr=np.sqrt(bkg_spectrum), ds=ds, label=f'bkg (#{bkg_doc["_id"]}, {bkg_start_date})')
            axs[0,1].errorbar(range(32), sig_spectrum  , yerr=np.sqrt(sig_spectrum), ds=ds, label='before bkg. sub.')
            spec_ebins=[seb.to_keV(i, i) for i in range(32)]

            axs[0,1].set_xlabel('Energy (keV)')
            axs[0,1].set_xticklabels(spec_ebins)
            axs[0,1].set_ylabel('Counts')
            axs[0,1].set_yscale('log')
            axs[0,1].set_xscale('log')
            axs[0,1].legend()


            cfl_tot_counts=np.sum(cfl_pixel_bkgsub_counts)
            pattern=(np.array(res['norm_counts'])*cfl_tot_counts).tolist()
            pattern_error=(np.array(res['norm_counts_errors'])*cfl_tot_counts).tolist()
            pattern.append(mean_fluence)
            pattern_error.append(mean_fluence_error)

            axs[1,0].errorbar(range(13), pattern, yerr=pattern_error,marker='o',
                    label="measurements")
            best_fit_pattern=(np.array(res['pixel_area_norm'])*cfl_tot_counts).tolist()
            best_fit_pattern.append(res['cfl_fluence'])
            best_fit_pattern_error=[0]*12
            best_fit_pattern_error.append(res['cfl_fluence_error'])
            axs[1,0].plot(range(13), best_fit_pattern, marker='o',
                    label='best match')
            axs[1,0].set_xlabel('Pixel ID')
            axs[1,0].set_ylabel('Counts')
            axs[1,0].set_title('CFL pattern')
            axs[1,0].legend()


            chi2 = res['delta_chi2'][1]
            x = np.linspace(res['x'][0] * 3600, res['x'][1] * 3600, res['x'][2])
            y = np.linspace(res['y'][0] * 3600, res['y'][1] * 3600, res['y'][2])
            X,Y=np.meshgrid(x,y)
            chi2_limit=22
            chi2_map[chi2_map>chi2_limit]=-1
            im=axs[1,1].pcolormesh(X,Y, np.transpose(chi2_map), cmap='RdPu', vmin=0)
            #fig.colorbar(im, ax=axs[1,1])
            axs[1,1].plot([xloc],[yloc],marker='+',markersize=10)
            axs[1,1].text(xloc+80,yloc-40, f'[{xloc:.1f}, {yloc:.1f}]', family='sans-serif', size=14)

            axs[1,1].grid()
            axs[1,1].set_xlim(-1600,1600)
            axs[1,1].set_ylim(-1600,1600)
            axs[1,1].set_xlabel('X (arcsec)')
            axs[1,1].set_ylabel('Y (arcsec)')
            axs[1,1].set_title('Flare location (CFL solution)')
            circle = Circle((0, 0), res['sun_angular_diameter']*0.5, alpha=0.8, ec='red', fc='none')
            axs[1,1].add_patch(circle)
            axs[1,1].set_aspect('equal')
            axs[1,1].text(0,res['sun_angular_diameter']*0.5,
                    f'North', family='sans-serif', size=12)
            plt.suptitle(f'STIX flare #{flare_id} (BSD #{_id})')


            
            if plot_fig:
                plt.show()
            if save_fig:
                try:
                    plt.savefig(fig_filename, dpi=300)
                    print('Fig saved to ', fig_filename)
                except Exception as e:
                    print(e)
            plt.close()
            """

    def process_L1_BSD_in_file(self, file_id):
        print("Processing file:", file_id)
        bsd_cursor = self.bsd_db.find({
            'run_id': file_id,
            'SPID': 54115
        }).sort('_id', 1)
        for doc in bsd_cursor:
            self.process_one_flare(doc, plot_fig=False, save_fig=True)

    def process_L1_BSD(self, bsd_id, bkg_bsd_id=None):
        doc = self.bsd_db.find_one({'_id': bsd_id, 'SPID': 54115})
        self.process_one_flare(doc, bkg_bsd_id)


if __name__ == '__main__':
    flp = FlareDataAnalyzer()
    if len(sys.argv) < 2:
        print('flare_location_solver run_id')
        print('flare_location_solver run_id_start id_end')

        flp.process_L1_BSD(2720, 2722)
    elif len(sys.argv) == 2:
        flp.process_L1_BSD_in_file(int(sys.argv[1]))
    else:
        for i in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
            flp.process_L1_BSD_in_file(i)
