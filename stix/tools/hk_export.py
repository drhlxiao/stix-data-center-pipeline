import pymongo
import sys
sys.path.append('.')
from stix.core import stix_datatypes as sdt
from stix.core import stix_datetime
connect = pymongo.MongoClient('localhost', 27017)
stix=connect['stix']
db=stix['packets']

HK={'NIXD0023': 'Instrument mode', 'NIXD0025': 'HK_DPU_PCB_T', 'NIXD0026': 'HK_DPU_FPGA_T', 'NIXD0027': 'HK_DPU_3V3_C', 'NIXD0028': 'HK_DPU_2V5_C', 'NIXD0029': 'HK_DPU_1V5_C', 'NIXD0030': 'HK_DPU_SPW_C', 'NIXD0031': 'HK_DPU_SPW0_V', 'NIXD0032': 'HK_DPU_SPW1_V', 'NIXD0038': 'HK_ASP_REF_2V5A_V', 'NIXD0039': 'HK_ASP_REF_2V5B_V', 'NIXD0040': 'HK_ASP_TIM01_T', 'NIXD0041': 'HK_ASP_TIM02_T', 'NIXD0042': 'HK_ASP_TIM03_T', 'NIXD0043': 'HK_ASP_TIM04_T', 'NIXD0044': 'HK_ASP_TIM05_T', 'NIXD0045': 'HK_ASP_TIM06_T', 'NIXD0046': 'HK_ASP_TIM07_T', 'NIXD0047': 'HK_ASP_TIM08_T', 'NIXD0048': 'HK_ASP_VSENSA_V', 'NIXD0049': 'HK_ASP_VSENSB_V', 'NIXD0050': 'HK_ATT_V', 'NIXD0051': 'HK_ATT_T', 'NIXD0052': 'HK_HV_01_16_V', 'NIXD0053': 'HK_HV_17_32_V', 'NIXD0054': 'DET_Q1_T', 'NIXD0055': 'DET_Q2_T', 'NIXD0056': 'DET_Q3_T', 'NIXD0057': 'DET_Q4_T', 'NIXD0035': 'HK_DPU_1V5_V', 'NIXD0036': 'HK_REF_2V5_V', 'NIXD0037': 'HK_DPU_2V9_V', 'NIXD0024': 'HK_PSU_TEMP_T', 'NIXD0001': 'SW Version Number', 'NIXD0002': 'CPU load', 'NIXD0003': 'Archive Memory usage', 'NIXD0166': 'Autonomous ASW boot stat', 'NIXD0167': 'Memory load ena flag', 'NIXD0004': 'IDPU identifier', 'NIXD0005': 'Active SPW link', 'NIXD0168': 'Overruns for tasks', 'NIXD0169': 'Watchdog state', 'NIXD0079': 'Received SpW packets', 'NIXD0078': 'Rejected SpW packets', 'NIXD0070': 'En/Dis Detector Status', 'NIXD0080': 'SPW1 - power status', 'NIXD0081': 'SPW0 - power status', 'NIXD0082': 'Q4 - power status', 'NIXD0083': 'Q3 - power status', 'NIXD0084': 'Q2 - power status', 'NIXD0085': 'Q1 - power status', 'NIXD0086': 'ASPECT B - power status', 'NIXD0087': 'ASPECT A - power status', 'NIXD0088': 'ATT M2 - moving', 'NIXD0089': 'ATT M1 - moving', 'NIXD0090': 'HV17-32 - enabled status', 'NIXD0091': 'HV01-16 - enabled status', 'NIXD0092': 'LV - enabled status', 'NIXD0066': 'HV1 depolar in progress', 'NIXD0067': 'HV2 depolar in progress', 'NIXD0068': 'ATT AB flag - OPEN', 'NIXD0069': 'ATT BC flag - CLOSED', 'NIX00072': 'Med value of trig acc', 'NIX00073': 'Max value of trig acc', 'NIXD0074': 'HV regulators mask', 'NIXD0077': 'TC(20,128) seq cnt', 'NIX00076': 'Attenuator motions', 'NIX00078': 'HK_ASP_PHOTOA0_V', 'NIX00079': 'HK_ASP_PHOTOA1_V', 'NIX00080': 'HK_ASP_PHOTOB0_V', 'NIX00081': 'HK_ASP_PHOTOB1_V', 'NIX00094': 'Attenuator currents', 'NIXD0075': 'HK_ATT_C', 'NIXD0058': 'HK_DET_C', 'NIX00085': 'FDIR function status'}

def request_hk_packets(start_utc, end_utc):
    names=HK.keys()
    start_unix=stix_datetime.utc2unix(start_utc)
    end_unix= stix_datetime.utc2unix(end_utc)
    cur=db.find({'header.SPID':54102, 'header.unix_time':{'$gt':start_unix,'$lt':end_unix}}).sort('header.unix_time',1)
    print(f'Number of packets:{cur.count()}')
    last_time=0
    csv_filename=f'stix_hk_{start_utc}_{end_utc}.csv'
    csv_filename=csv_filename.replace(':','')
    fcsv=open(csv_filename,'w')
    header_name='Packet ID, UTC,'+','.join(HK.values())+'\n'
    fcsv.write(header_name)
    print('request data...')
    i=0
    indexes=[]
    for pkt in cur:
        header=pkt['header']
        parameters=pkt['parameters']
        if header['unix_time']==last_time:
            continue
        if i==0:
            for name in names:
                for j, p in enumerate(parameters):
                    if p[0]==name:
                        indexes.append(j)
                        break



        last_time=header['unix_time']
        _id=pkt['_id']
        line=f'{_id},{header["UTC"]},'
        line+=','.join([str(parameters[ix][2]) for ix in indexes])
        line+='\n'
        fcsv.write(line)
        i+=1
        if i%1000==0:
            print(i)

    print('Done.')

    

request_hk_packets('2020-04-15T00:00:00','2021-05-04T00:00:00')







        
        
        







