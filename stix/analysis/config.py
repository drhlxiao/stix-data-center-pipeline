#!/usr/bin/python
'''
get stix configuration from the database

'''
import sys
sys.path.append('.')
import numpy as np
from stix.core import mongo_db as sdb
from stix.spice import datetime as sdt
mdb = sdb.MongoDB()
scdb=mdb.get_collection('config')
caldb=mdb.get_collection('calibration_runs')
MIN_CALIBRATION_DURATION=12*3600
NOMINAL_EBIN_EDGES= np.array([0,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36,
    40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150, np.inf
    ])
#NOMINAL_EBINS=[(ELUT_EBINS_ENERGIES[i],ELUT_EBINS_ENERGIES[i+1]) for i in range(32)]


class Elut(object):
    def __init__(self,utc):
        self.utc=utc
        self.onboard_elut=Elut.get_onboard_elut(utc)
        self.calibration_run_elut=Elut.get_calibration_run_elut(utc)

        #adc = slope * E + offset 
        oslp,oofs=np.array(self.onboard_elut['slopes']), np.array(self.onboard_elut['offsets'])
        cslp,cofs=np.array(self.calibration_run_elut['slopes']), np.array(self.calibration_run_elut['offsets'])
        f=lambda energy: (oslp*energy+oofs-cofs)/cslp  # true energies  for 384 pixels, 
        #self.true_energy_bin_edges=np.array([ f(x) for x in NOMINAL_EBIN_EDGES]) # dimension 384 x 33 
        self.true_energy_bin_edges=np.round(np.array([ f(x) for x in NOMINAL_EBIN_EDGES]),4) # dimension 384 x 33 
    def get_data(self):
        return {'data':{
                    'onboard':self.onboard_elut, 
                    'calibration':self.calibration_run_elut,
                    'true_energy_bin_edges':self.true_energy_bin_edges,
                    'energy_bin_edges':NOMINAL_EBIN_EDGES,
                },
                'info': self.info(),
                }
    def info(self):
        return {
                'calibration_run':{
                    'obs_begin': self.calibration_run_elut['start_utc'],
                    'duration': self.calibration_run_elut['duration'],
                    'run_id': self.calibration_run_elut['run_id']
                    },
                'onboard_elut':{
                    'upload_time_range':self.onboard_elut['upload_time_range'],
                    }
                }
    def get_calibration_data(self):
        return self.calibration_run_elut
    #def get_onboard_energy_lut(self):
    #    return self.onboard_elut
    def get_true_ebin_edges(self):
        return self.true_energy_bin_edges
    def get_pixel_true_ebins(self, pixel):
        ebins=self.true_energy_bin_edges[:,pixel] #retrieve the column  
        return np.array([ (ebins[i], ebins[i+1]) for i in range(32)])

    @staticmethod
    def get_onboard_elut(utc):
        unix=sdt.utc2unix(utc)
        elut={}
        min_time=5e9
        max_time=0
        #pkt_ids=[]
        offsets=[0]*384
        slopes=[0]*384
        for i in range(384):
            pixel_elut=list(scdb.find({'pixel_id':i, 'type':'elut', 'execution_unix':{'$lte':unix}}).sort('execution_unix',-1).limit(1))
            if pixel_elut:
                offsets[i]=pixel_elut[0]['offset']
                slopes[i]=pixel_elut[0]['slope']
                uptime=pixel_elut[0]['execution_unix']
                if uptime<min_time:
                    min_time=uptime
                if uptime>max_time:
                    max_time=uptime
            #pkt_ids.append(pixel_elut[0]['packet_id'])
        elut={'slopes':slopes,
                'offsets':offsets,
                'upload_time_range':[sdt.unix2utc(min_time), sdt.unix2utc(max_time)],
                'energy_bin_edges':NOMINAL_EBIN_EDGES,
                #'packet_ids':pkt_ids
                }
        return elut
    @staticmethod
    def get_calibration_run_elut(utc):
        unix=sdt.utc2unix(utc)
        run=list(caldb.find({'start_unix_time':{'$lte':unix}, 'analysis_report':{'$exists':True},
            'duration':{'$gt':MIN_CALIBRATION_DURATION}}).sort('start_unix_time',-1).limit(1))
        res={}
        if run:
            res={
                    'slopes':np.round(run[0]['analysis_report']['slope'],4),
                    'offsets':np.round(run[0]['analysis_report']['offset'],4),
                    'slope_errors':np.round(run[0]['analysis_report']['slope_error'],4),
                    'offset_errors':np.round(run[0]['analysis_report']['offset_error'],4),
                    'run_id':run[0]['_id'],
                    'duration':run[0]['duration'],
                    'start_unix_time':run[0]['start_unix_time'],
                    'start_utc':sdt.unix2utc(run[0]['start_unix_time'])
                    }
        return res




if __name__=='__main__':
    print(Elut.get_onboard_elut('2021-08-31T00:00:00'))
    print(Elut.get_calibration_run_elut('2021-08-31T00:00:00'))
    t=Elut('2021-08-31T00:00:00')
    print('=======')
    print(t.info())
    print("Pixel 0 ELUT:")
    print(t.get_pixel_true_ebins(0))
    print("Pixel 12 ELUT:")
    print(t.get_pixel_true_ebins(12))









