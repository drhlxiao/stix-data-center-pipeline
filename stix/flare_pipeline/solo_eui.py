#Author Hualin Xiao (hualin.xiao@fhnw.ch)
# Download data from Solar Orbiter Data Archive
import requests 
import wget
import tempfile
from dateutil import parser as dtparser
from datetime import timedelta
from datetime import datetime
from tqdm import tqdm
from sunpy.map import Map
from pprint import pprint 
import re

prefix='https://wwwbis.sidc.be/EUI/data/L1'
esa_url='http://soar.esac.esa.int/soar-sl-tap/tap/sync'        

def utc2datetime(utc):
    if not utc.endswith('Z'):
        utc += 'Z'
    return dtparser.parse(utc)

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
            fdt=utc2datetime(f)
            fits_time.append([fdt,url, fname])
    return fits_time
def get_eui_fits_list(start, end):
    sdt=utc2datetime(start)
    edt=utc2datetime(end)
    num=(edt-sdt).days+1
    results=[]
    for i in range(num):
        dt=sdt+timedelta(days=i)
        fits=find_eui_data(dt)
        for entry in fits:
            if entry[0]>=sdt and entry[0]<=edt:
                results.append(entry)
    return results


def query_esa_tap_server(start_utc, end_utc, instrument='EUI', level='L2'):
    sdt=utc2datetime(start_utc)
    edt=utc2datetime(end_utc)
    start_utc=sdt.strftime('%Y-%m-%d %H:%M:%S.000')
    end_utc=edt.strftime('%Y-%m-%d %H:%M:%S.000')

    sql=f'''SELECT * FROM v_sc_data_item WHERE data_type='SCI'  AND   (end_time >'{start_utc}')   AND   (begin_time <'{end_utc}')   AND   ( (instrument ='EUI') )   AND   ( (level ='L2') )'''

    data={'REQUEST': 'doQuery',
            'LANG': 'ADQL',
            'FORMAT': 'JSON',
            'PHASE': 'RUN',
            'QUERY': sql, 
            'PAGE': 1,
            'PAGE_SIZE': 500}
    print(sql)
    x = requests.post(esa_url, data = data)
    return x.json()

def get_url(row):
    file_id=row[1]
    return f'http://soar.esac.esa.int/soar-sl-tap/data?retrieval_type=PRODUCT&data_item_id={file_id}&product_type=SCIENCE'
def plot_eui(folder, start_utc, end_utc, peak_utc, instrument='EUI'):
    resp=query_esa_tap_server(start_utc, end_utc, instrument=instrument)
    pprint(resp)
    if not resp['data']:
        print("eui not available")
        return None
    peak_ut=utc2datetime(peak_utc).timestamp()
    eui_fsi174 = { abs(utc2datetime(row[0]).timestamp() - peak_ut):row  for row in resp['data'] if 'FSI174' in row[4]}
    eui_fsi304= { abs(utc2datetime(row[0]).timestamp() - peak_ut):row  for row in resp['data'] if 'FSI304' in row[4]}
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

    return  map_fsi304, map_fsi174


if __name__=='__main__':
    plot_eui('.', '2022-03-01T00:00:00','2022-03-01T02:10:00','2022-03-01T00:05:00',
            'EUI')
