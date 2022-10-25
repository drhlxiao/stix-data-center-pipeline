#Author Hualin Xiao (hualin.xiao@fhnw.ch)
# Download data from Solar Orbiter Data Archive
import re
import sys
import os
import wget
import requests 
import tempfile
from dateutil import parser as dtparser
from datetime import timedelta
from datetime import datetime
from stix.spice import time_utils as stu
from tqdm import tqdm
from sunpy.map import Map
from pprint import pprint 
import sunpy
from matplotlib import pyplot as plt
import matplotlib
from sunpy.map import make_fitswcs_header
from sunpy.coordinates.frames import HeliocentricEarthEcliptic, HeliographicStonyhurst
from sunpy.visualization import colormaps as scm
from stix.core import logger
from stix.core import mongo_db as db
from matplotlib import cm
matplotlib.use('Agg')

logger = logger.get_logger()
mdb = db.MongoDB()
flare_image_db = mdb.get_collection('flare_images')
cmap_fsi174='Greys'

prefix='https://wwwbis.sidc.be/EUI/data/L1'
esa_url='http://soar.esac.esa.int/soar-sl-tap/tap/sync'        

def find_eui_data(dt):
    url=f'{prefix}/{dt.year}/{dt.month:02d}/{dt.day:02d}'
    html=requests.get(url).text
    if not html:
        return []
    fits= re.findall(r'href=\"(solo_L1_eui.*fits)"', html)
    fits_time=[]
    for fname in fits:
        ft= re.findall(r'\d{8}T\d{6}', fname)
        for f in ft:
            fdt=stu.utc2datetime(f)
            fits_time.append([fdt,url, fname])
    return fits_time
def get_eui_fits_list(start, end):
    sdt=stu.utc2datetime(start)
    edt=stu.utc2datetime(end)
    num=(edt-sdt).days+1
    results=[]
    for i in range(num):
        dt=sdt+timedelta(days=i)
        fits=find_eui_data(dt)
        for entry in fits:
            if entry[0]>=sdt and entry[0]<=edt:
                results.append(entry)
    return results


def query_esa_tap_server(start_dt, end_dt, instrument='EUI', level='L2'):
    start_utc=start_dt.strftime('%Y-%m-%d %H:%M:%S.000')
    end_utc=end_dt.strftime('%Y-%m-%d %H:%M:%S.000')

    sql=f'''SELECT * FROM v_sc_data_item WHERE data_type='SCI'  AND   (end_time >'{start_utc}')   AND   (begin_time <'{end_utc}')   AND   ( (instrument ='EUI') )   AND   ( (level ='L2') )'''

    data={'REQUEST': 'doQuery',
            'LANG': 'ADQL',
            'FORMAT': 'JSON',
            'PHASE': 'RUN',
            'QUERY': sql, 
            'PAGE': 1,
            'PAGE_SIZE': 500}
    x = requests.post(esa_url, data = data)
    return x.json()

def get_url(row):
    file_id=row[1]
    return f'http://soar.esac.esa.int/soar-sl-tap/data?retrieval_type=PRODUCT&data_item_id={file_id}&product_type=SCIENCE'
def get_eui_maps(start_utc, end_utc, peak_utc=None, instrument='EUI', padding = 300):

    res=[]
    start_dt=stu.utc2datetime(start_utc) - timedelta(seconds=padding)
    end_dt=stu.utc2datetime(end_utc) + timedelta(seconds=padding)


    resp=query_esa_tap_server(start_dt, end_dt, instrument=instrument)
    desc=[]
    if not resp['data']:
        print("eui not available")
        return [], []
    if peak_utc:
        peak_ut=stu.utc2datetime(peak_utc).timestamp()
    else:
        peak_ut=(start_dt.timestamp() + end_dt.timestamp())/2.

    eui_fsi174 = { abs(stu.utc2datetime(row[0]).timestamp() - peak_ut):row  for row in resp['data'] if 'FSI174' in row[4]}
    eui_fsi304= { abs(stu.utc2datetime(row[0]).timestamp() - peak_ut):row  for row in resp['data'] if 'FSI304' in row[4]}
    temp_folder=tempfile.gettempdir()
    map_fsi304, map_fsi174=None, None
    if eui_fsi174:
        url_fsi174=get_url(eui_fsi174[min(eui_fsi174.keys())])
        filename_fsi174=wget.download(url_fsi174, out=temp_folder)
        map_fsi174=Map(filename_fsi174)
        os.unlink(filename_fsi174)
    if eui_fsi304:
        url_fsi304=get_url(eui_fsi304[min(eui_fsi304.keys())])
        filename_fsi304=wget.download(url_fsi304, out=temp_folder)
        map_fsi304=Map(filename_fsi304)
        os.unlink(filename_fsi304)
    if map_fsi174: 
        res.append(map_fsi174)
        desc.append('FSI174')
    if map_fsi304: 
        res.append(map_fsi304)
        desc.append('FSI304')

    return  res, desc

def process_one(image_run_id):
    logger.info(f'Processing: {image_run_id}')
    doc=flare_image_db.find_one({'_id':image_run_id})
    process_one_doc(doc)

def process_between(start, end):
    docs=flare_image_db.find({'_id':{'$gt':start, '$lt':end} } ).sort('_id',1)
    for doc in docs:
        logger.info(f'processing STIX image {doc["_id"]}')
        try:
            process_one_doc(doc)
        except Exception as e:
            logger.error(f'Failed when processing STIX image {doc["_id"]}')
            logger.error(str(e))
def process_one_doc(doc):
    if not doc:
        logger.info(f'Could not create EUI image, no flare image {doc["_id"]}')
        return

    try:
        image_clean=doc['fits']['image_clean']
    except (KeyError, TypeError):
        logger.info('Could not create EUI image, can not read info from STIX EM image !')
        return
   
    def exists(objs, ts):
        for o in objs:
            if o['type'] == ts:
                return True
        return False

    stix_map = sunpy.map.Map(image_clean)

    stix_bottom_left = stix_map.bottom_left_coord
    stix_top_right = stix_map.top_right_coord

    report=doc.get('report',{})
    out_folder=doc['idl_config']['folder']
    fout_prefix=doc["idl_config"]["prefix"]

    start_utc=stu.unix2utc(doc['start_unix'])
    end_utc=stu.unix2utc(doc['end_unix'])

    try:
        res, desc = get_eui_maps(start_utc, end_utc)
    except Exception as e:
        logger.error(str(e))
        return
    if not res:
        logger.info('No EUI observations')
        return


    report=doc.get('report',{})
    gcolor='w'


    for m, wavlen  in zip(res, desc):
        #m.peek()
        fig = plt.figure(figsize=(5, 5))
        msub = m.submap(stix_bottom_left.transform_to(m.coordinate_frame), top_right=stix_top_right.transform_to(m.coordinate_frame))
        comp_map=sunpy.map.Map( msub, stix_map, composite=True)
        ps = comp_map.get_plot_settings(0)
        ps['cmap'] = plt.get_cmap('sdoaia171').reversed() if '174' in wavlen else cm.hot
        comp_map.set_plot_settings(0, ps)
        levels = [50, 60, 70,80,90] #contour levels in percent
        comp_map.set_levels(index=1, levels=levels, percent=True)
        ax = fig.add_subplot(111, projection=stix_map)
        erange=f"{doc['energy_range'][0]} - {doc['energy_range'][1]} keV "
        title=f"EUI {wavlen} {m.meta['date-obs'] } \n STIX CLEAN {erange} {stix_map.meta['date_avg']}  "
        comp_map.plot(axes=ax, title=title)

        fname=os.path.join(out_folder, f'{fout_prefix}_eui_{wavlen}.png')
        eui_croped_fits=os.path.join(out_folder, f'{fout_prefix}_eui_{wavlen}.fits')
        fig.savefig(fname)
        msub.save(eui_croped_fits, overwrite=True)
        report[f'V_eui-{wavlen}']={'filename':fname,
                'fits':eui_croped_fits,
                'title': f'Solar Orbiter EUI {wavlen} image'}
        

    updates={'report':report}
    flare_image_db.update_one({'_id': doc['_id']},{'$set':updates})



if __name__=='__main__':
    if len(sys.argv)==2:
        process_one(int(sys.argv[1]))
    elif len(sys.argv)==3:
        process_between(int(sys.argv[1]),int(sys.argv[2]))
    else:
        print('process <id> or process start  end_id')
