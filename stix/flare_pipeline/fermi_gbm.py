import sys
import math
import numpy as np
import wget
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
from astropy.io  import fits 
from datetime import datetime, timedelta
from gbm.finder import TriggerFtp, ContinuousFtp
from gbm.data import TTE,Ctime, RSP, PosHist
from gbm import time as gt
from gbm.binning.unbinned import bin_by_time
from gbm.detectors import Detector
from gbm import coords as gc
import matplotlib.pyplot as plt
from gbm.plot import Lightcurve, Spectrum
import dateutil.parser as dtp

from gbm.background import BackgroundFitter
from gbm.background.binned import Polynomial

from stix.core import mongo_db as db
from stix.spice.stix_datetime import unix2utc
mdb = db.MongoDB()
db=mdb.get_collection('flares')

#from astropy.utils.data import download_file
#https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs/notebooks/TteData.html
#https://fermi.gsfc.nasa.gov/ssc/library/support/Science_DP_ICD_RevA.pdf
#API https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs/api/api-data.html#gbm.data.Ctime

def to_met(t):
    if not t:
        return None
    dt=None
    if isinstance(t,str):
        dt=dtp.parse(t)
    elif isinstance(t,datetime):
        dt=t
    
    if dt is not None:
        return gt.Met.from_datetime(dt).met.astype(float)
    return None
def to_datetime(t):
    if not t:
        return None
    dt=None
    if isinstance(t,str):
        dt=dtp.parse(t)
    elif isinstance(t,datetime):
        dt=t
    return dt
        
def download_fermi(dt: datetime, folder='.', file_type='poshist', detector=''):
    d=dt.date()
    hour=dt.hour
    if file_type=='poshist':
        fname=f'glg_poshist_all_{d.year-2000:02d}{d.month:02d}{d.day:02d}_v00.fit'
    elif file_type=='tte' and detector:
        fname=f'glg_tte_{detector}_{d.year-2000:02d}{d.month:02d}{d.day:02d}_{hour:02d}z_v00.fit.gz'
    else:
        return None
    local_file=Path(folder)/fname
    if not local_file.is_file() and fname:
        
        url= f'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/{d.year}/{d.month:02d}/{d.day:02d}/current/{fname}'
        print('Downloading ', url)
        try:
            local_filename=wget.download( url,folder)
        except:
            print('Failed to download: ', url)
            return None
    if local_file.is_file():
        return local_file
    return None
    
def download_fermi_sunward_detector_data( start_utc, end_utc,folder='.', download_tte=True):
    """ 
        
        download fermi  TTE and poshis files for sunward detectors in the given time frame
        for light curve and spectral analysis
        Arguments
        
        start_utc: str or datetime.datetime
            start time 
        end_utc:
            str or datetime.datetime
        Returns:
            res:  dict or None
                contains time detector name and downloaded tte filename
                Return None if not found
    """
    detector_names = ('n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb')#, 'b0', 'b1')
    start_dt=to_datetime(start_utc)
    end_dt=to_datetime(end_utc)
    dates=set([start_dt.date(), end_dt.date()])
    #same_day=True if len(dates)==1 else False
    duration=(end_dt-start_dt).total_seconds()
    dt_step=int(duration/600)
    
    delta_times=[timedelta(seconds=sec) for sec in np.arange(0,duration,dt_step)]
    last_date=None
    result={'time':[],'detector':[],'angle':[],'tte':[]}
    ph=None
    last_detector=None
    for delta_t in delta_times:
        dt=start_dt+delta_t
        
        if dt.date()!=last_date:
            local_file=download_fermi(dt, folder, file_type='poshist', detector='')
            if not local_file:
                return result
            ph=PosHist.open(local_file)
        last_date=dt.date()
        mt=gt.Met.from_datetime(dt).met        
        sun_vis=ph.get_sun_visibility(mt.astype(float))
        if sun_vis:
            sun_ra, sun_dec=gc.get_sun_loc(mt)
            min_angle=math.inf
            det_sunward=None
            for det in detector_names:
                angle=ph.detector_angle(sun_ra, sun_dec, det, mt)
                if angle<min_angle:
                    min_angle=angle
                    det_sunward=det
            if download_tte:    
                tte_file=  download_fermi(dt, folder, 'tte', det_sunward)
                result['tte'].append(tte_file)

            result['time'].append(gt.Met(mt).datetime) 
            result['detector'].append(det_sunward)
            result['angle'].append( min_angle)
            
            
    return result
def plot_fermi_gdm_lightcurve_spectrum(start_utc='', end_utc='', bkg_start=None, bkg_end=None, energy_range=(10,26), image_folder='.', cache_folder='.', flare_id='', entry_id=''):
    res=download_fermi_sunward_detector_data( start_utc, end_utc,cache_folder)
    #print(res)
    if not res['tte']:
        print('No sun pointing detector  avaible')
    tte_flist=set(res['tte'])
    detectors=set(res['detector'])
    tte_list=[]
    for fname in tte_flist:
        print('Opening ', fname)
        t = TTE.open(fname)
        tte_list.append(t)
    if len(tte_list)>1:
        tte=TTE.merge(tte_list)
    elif len(tte_flist)==1:
        tte = tte_list[0]
    else:
        return None
    
    print(tte.time_range)
    print(tte.energy_range)
    start_met=to_met(start_utc)
    end_met=to_met(end_utc)
    if not start_met or not end_met:
        return None
    time_sliced_tte = tte.slice_time([start_met, end_met])
    energy_sliced_tte = time_sliced_tte.slice_energy(energy_range)
    phaii =energy_sliced_tte.to_phaii(bin_by_time, 4, time_ref=start_met)
    fig, axs=plt.subplots()
    lcplot = Lightcurve(data=phaii.to_lightcurve(), axis=axs)
    #plt.show()
    #spectrum = time_sliced_tte.to_spectrum()
    #specplot = Spectrum(data=spectrum, axis=axs[1])
    det='-'.join(detectors)
    axs.set_title(f'GBM det {det} obs for flare {flare_id} ({entry_id})')
    axs.set_ylabel('Start at {start_utc}')
    

    #plt.figtext(0.5, 0.01, f'{detectors}  {start_utc} -- {end_utc} E: {energy_range} keV', 
    #            ha="left", )
    fname=Path(image_folder)/f'fermi_{entry_id}_{det}_{flare_id}.png'
    fig.savefig(fname)
    return fname
    
def plot_fermi_for_flare(image_folder, cache_folder,_id, overwrite=False):
    #_id is the flare db entry number
    print('Processing flare #', _id)
    flare_doc=db.find_one({'_id':_id})
    if not flare_doc:
        print(f'Flare {_id} does not exist!')
        return None
    flare_id=flare_doc['flare_id']
    if flare_doc['peak_counts']<600:
        print(f'Flare #{flare_id} skipped because its counts too low ')
        return None
        

    key='fermi'
    if  mdb.get_flare_pipeline_products(_id, key) and overwrite == False:
        print(f'Fermi for Flare #{_id} has been created!')
        return None
    if flare_doc['goes']['class']=='A':
        print(f'Flare #{flare_id} skipped because it was not seen by GOES')
        return None
    light_time_diff=flare_doc['ephemeris']['light_time_diff']
    start=flare_doc['start_unix']+light_time_diff-120
    end=flare_doc['end_unix'] +light_time_diff + 120
    start_utc=unix2utc(start)
    end_utc=unix2utc(end)
    fname=plot_fermi_gdm_lightcurve_spectrum(start_utc, end_utc, image_folder=image_folder, cache_folder=cache_folder, flare_id=flare_id, 
            entry_id=_id)
    if fname:
        mdb.update_flare_pipeline_products(_id, key, [str(fname)])
    #plot_fermi_gdm_lightcurve_spectrum('2021-10-07T02:40:00', '2021-10-07T02:50:00')

if __name__=='__main__':
    image_folder='/data/fermi/images'
    cache_folder='/data/fermi/fits'
    ids=[]
    if len(sys.argv)==1:
        print('fermi_gbm <flare_entry_id> <end_id>')
        plot_fermi_for_flare('/tmp','/tmp',2303)

    elif len(sys.argv)==2:
        ids.append(int(sys.argv[1]))
    elif len(sys.argv)==3:
        ids=range(int(sys.argv[1]), int(sys.argv[2]))

    for _id in ids:
        try:
            plot_fermi_for_flare(image_folder, cache_folder, _id)
        except Exception as e:
            print(e)


    
