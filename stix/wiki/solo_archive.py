import requests 
from dateutil import parser as dtparser
from datetime import timedelta
from datetime import datetime
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



#request data directly for esa archive
def query_esa_tap_server(start_utc, end_utc, instrument='EUI', level='L1'):
    sdt=utc2datetime(start_utc)
    edt=utc2datetime(end_utc)
    start_utc=sdt.strftime('%Y-%m-%d %H:%M:%S.000')
    end_utc=edt.strftime('%Y-%m-%d %H:%M:%S.000')
    data={'REQUEST': 'doQuery',
            'LANG': 'ADQL',
            'FORMAT': 'JSON',
            'PHASE': 'RUN',
            'QUERY': f"SELECT data_item_oid, data_item_id,  descriptor,  begin_time, end_time,  filename, filesize FROM v_sc_data_item WHERE data_type='SCI'  AND   (begin_time >'{start_utc}')   AND   (begin_time <'{end_utc}')   AND   ( (instrument ='{instrument}') )   AND   ( (level ='{level}') )  ORDER BY begin_time ASC", 
            'PAGE': 1,
            'PAGE_SIZE': 500}
    x = requests.post(esa_url, data = data)
    return x.json()
def request_solo_payload_data(instrument):
    data=[]
    levels={'eui':'L1',
                'epd':'L2',
                }
    payloads={'eui':'EUI',
                'epd':'EPD'}
    data=tap.query_esa_tap_server(start, end, payloads[instrument], levels[instrument])
    return data


