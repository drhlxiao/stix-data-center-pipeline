import os
import sys
import numpy as np
from scipy import interpolate
from stix.spice import time_utils
from stix.core import datatypes as sdt
from stix.core.energy_bins import StixEnergyBins
from stix.core import mongo_db as db
from matplotlib import pyplot as plt

QLLC_SPID = 54118
QLLC_BKG_SPID = 54119
QLBKG_SPID= 54120

mdb = db.MongoDB()

def set_db(db):
    global mdb
    mdb=db

def extrapolate_counts(lc_times, lc_counts, tx):
    max_index=np.argmax(lc_counts)
    lc_times=lc_times[:max_index]
    lc_counts=lc_counts[:max_index]
    #make sure we use the counts before the max counts

    f = interpolate.interp1d(lc_times,lc_counts, fill_value='extrapolate')
    return f(tx)


class QLAnalyzer(object):
    def __init__(self):
        self._data=None

    @property
    def data(self):
        return self._data

    def get_counts_range(self, data=None):
        if data is None and self._data is not None:
            data=self._data
        if not data:
            print('Data not ready')
            return {}

        num_lcs=len(data['lcs'].keys())
        results=np.zeros((2,num_lcs))
        for ilc in range(num_lcs):
            counts=data['lcs'][ilc]
            results[0][ilc], results[1][ilc]=(np.min(counts),np.max(counts))
        return  results
    def peek(self, ax=None, ilcs=[0,1,2,3,4]):
        if self._data is None:
            print('Data not ready')
            return
        data=self._data
        if ax is None:
            fig, ax=plt.subplots()

        ql_times=time_utils.unix2datetime(data['time'])
        for i in ilcs:
            ax.plot(ql_times, data['lcs'][i], label=f'QL LC {data["energy_bins"]["names"][i]}')
        if 'rcr' in data:
            ax.plot(ql_times, data['rcr'], label=f'RCR')
        ax.legend()
        ax.set_yscale('log')
        return ax       


class LightCurveAnalyzer(QLAnalyzer):

    def __init__(self, data):
        self._data=data
        self._model=None
        self._model_fname=None

    #def set_att_in_pred_model(self,filename):
    #    """
    #    set random forest prediction model to predict counts when ATT is inserted
    #    """
    #    import joblib
    ##    self._model_fname=filename
    #    self._model=joblib.load(filename)


    @classmethod
    def from_database(cls, start_utc, end_utc):
        start_unix, end_unix=time_utils.anytime(start_utc,'unix'), time_utils.anytime(end_utc,'unix') 
        packets=mdb.get_quicklook_packets('lc', start_unix, end_unix-start_unix, sort_field='_id')
        data=LightCurveAnalyzer.parse(packets, True,start_unix, end_unix)
        return cls(data)
 

    @staticmethod
    def parse(packets, exclude_duplicated=True,start_unix=None, end_unix=None):
        #if not packets:
        #    return None
        lightcurves = {}
        unix_time = []
        energy_bins={}
        rcrs=[]
        triggers=[]
        last_time=0
        for pkt in packets:
            packet = sdt.Packet(pkt)
            if not packet.isa(QLLC_SPID):
                continue
            #fig = None

            scet_coarse = packet[1].raw
            scet_fine = packet[2].raw
            start_scet = scet_coarse + scet_fine / 65536.

            if start_scet<=last_time and exclude_duplicated:
                continue
            last_time=start_scet
            int_duration = (packet[3].raw + 1) * 0.1

            detector_mask = packet[4].raw
            pixel_mask = packet[6].raw

            num_lc = packet[17].raw

            compression_s = packet[8].raw
            compression_k = packet[9].raw
            compression_m = packet[10].raw

            num_lc_points = packet.get('NIX00270/NIX00271')[0]
            lc = packet.get('NIX00270/NIX00271/*.eng')[0]
            rcr = packet.get('NIX00275/*.raw')
            trig= packet.get('NIX00273/*.eng')[0]
            UTC = packet['header']['UTC']
            rcrs.extend(rcr[0])
            triggers.extend(trig)
            for i in range(len(lc)):
                if i not in lightcurves:
                    lightcurves[i] = []
                lightcurves[i].extend(lc[i])
            unix_time.extend([
                time_utils.scet2unix(start_scet + x * int_duration)
                for x in range(num_lc_points[0])
                ])

            if not energy_bins:
                energy_bin_mask= packet[16].raw
                tstart=time_utils.scet2unix(start_scet)
                energy_bins=StixEnergyBins.get_emask_energy_bins(energy_bin_mask, tstart)

        if not lightcurves:
            return None
        start_unix=unix_time[0] if start_unix is None else start_unix
        end_unix=unix_time[-1] if end_unix is None else end_unix

        unix_time=np.array(unix_time)
        sel=np.where( (unix_time>=start_unix) & (end_unix>=unix_time))

        utime=unix_time[sel]
        return {'time': utime,
                't':utime-utime[0],
                'lcs': {x: np.array(lightcurves[x])[sel] for x in lightcurves}, 
                'rcr':np.array(rcrs)[sel],
                'triggers':np.array(triggers)[sel],
                'energy_bins':energy_bins,'num':utime.size,'start_unix': utime[0],'end_unix':utime[-1]}

class BackgroundReportAnalyzer(QLAnalyzer):
    @classmethod
    def parse(cls, packets):
        if not packets:
            return None
        results=[]
        for pkt in packets:
            packet = sdt.Packet(pkt)
            if not packet.isa(QLBKG_SPID):
                continue
            scet_coarse = packet[1].raw
            scet_fine = packet[2].raw
            start_unix=time_utils.scet2unix(scet_coarse, scet_fine)
            tbin= (packet[3].raw + 1)*0.1
            num_samples=packet[14].raw
            samples=packet[14].children
            for i in range(num_samples):
                detector=samples[0+35*i][1]
                spectra=[samples[j+35*i][2] for j in range(32)]
                trig=samples[33+35*i][2]
                t=start_unix+samples[34+35*i][1]*tbin#number of integrations after the first one
                struct={
                        'detector':detector,
                        'spectra':spectra,
                        'triggers':trig,
                        'obs_time': t,
                        'tbin':tbin,
                        'packet_id':pkt['_id']
                        }
                results.push(struct)
        return results


class BKGLigthCurveAnalyzer(QLAnalyzer):

    def __init__(self, data):
        self._data=data

    @staticmethod
    def parse(packets, exclude_duplicated=True,start_unix=None, end_unix=None):
        #if not packets:
        #    return None
        lightcurves = {}
        unix_time = []
        energy_bins={}
        rcrs=[]
        triggers=[]
        last_time=0
        for pkt in packets:
            packet = sdt.Packet(pkt)
            if not packet.isa(QLLC_BKG_SPID):
                continue
            #fig = None

            scet_coarse = packet[1].raw
            scet_fine = packet[2].raw
            start_scet = scet_coarse + scet_fine / 65536.

            if start_scet<=last_time and exclude_duplicated:
                continue
            last_time=start_scet
            int_duration = (packet[3].raw + 1) * 0.1

            #detector_mask = packet[4].raw
            #pixel_mask = packet[6].raw

            num_lc = packet[14].raw

            compression_s = packet[11].raw
            compression_k = packet[12].raw
            compression_m = packet[13].raw

            num_lc_points = packet.get('NIX00273')[0]
            lc = packet.get('NIX00270/NIX00277/*.eng')[0]
            trig= packet.get('NIX00273/*.eng')[0]
            #rcr = packet.get('NIX00275/*.raw')
            UTC = packet['header']['UTC']
            triggers.extend(trig)
            for i in range(len(lc)):
                if i not in lightcurves:
                    lightcurves[i] = []
                lightcurves[i].extend(lc[i])
            unix_time.extend([
                time_utils.scet2unix(start_scet + x * int_duration)
                for x in range(num_lc_points)
                ])

            if not energy_bins:
                energy_bin_mask= packet[9].raw
                tstart=time_utils.scet2unix(start_scet)
                energy_bins=StixEnergyBins.get_emask_energy_bins(energy_bin_mask, tstart)

        if not lightcurves:
            return None
        unix_time=np.array(unix_time)

        start_unix=unix_time[0] if start_unix is None else start_unix
        end_unix=unix_time[-1] if end_unix is None else end_unix

        sel=np.where( (unix_time>=start_unix) & (end_unix>=unix_time))
        utime=unix_time[sel]
        return {'time': utime, 
                'time_since_T0':utime-utime[0],
                'lcs': {x: np.array(lightcurves[x])[sel] for x in lightcurves}, 'triggers':np.array(triggers)[sel],
                'energy_bins':energy_bins,'num':len(utime),'start_unix': utime[0],'end_unix':utime[-1]}



    @classmethod
    def from_database(cls, start_utc, end_utc):
        start_unix, end_unix=time_utils.anytime(start_utc,'unix'), time_utils.anytime(end_utc,'unix') 
        packets=mdb.get_quicklook_packets('bkg', start_unix, end_unix-start_unix, sort_field='_id')
        data=BKGLigthCurveAnalyzer.parse(packets, True,start_unix,end_unix)
        return cls(data)




class LightCurveMerger(QLAnalyzer):
    def __init__(self,qllc={}, qlbkg={}):
        self.ql_lc=qllc
        self.ql_bkg_lc=qlbkg
        self._data={'BKG_LC':qlbkg.data,
                'QL_LC': qllc.data
                }




    def peek(self, ax=None, ilcs=[0,1,2,3,4],):
        if self._data is None:
            print('Data not ready')
            return
        data=self._data
        #bkg_times=time_utils.unix2datetime(data['BKG_LC']['time'])
        if data['BKG_LC'] is None or data['QL_LC'] is None:
            print('Data is None')
            return
        if ax is None:
            fig, ax=plt.subplots()

        bkg_dt=time_utils.unix2datetime(data['BKG_LC']['time'])
        ql_dt=time_utils.unix2datetime(data['QL_LC']['time'])
        for i in ilcs:
            ax.plot(bkg_dt, data['BKG_LC']['lcs'][i], label=f'BKG LC {data["BKG_LC"]["energy_bins"]["names"][i]}')
            ax.plot(ql_dt, data['QL_LC']['lcs'][i], label=f'QL LC {data["QL_LC"]["energy_bins"]["names"][i]}')
        if 'rcr' in data['QL_LC']:
            ax.plot(ql_dt, data['QL_LC']['rcr'], label=f'RCR')
        ax.legend()
        ax.set_yscale('log')
        return ax
    def correct_att_in_counts(self,tgap=20):
        """
         here we should ensure the time selected should be between
         20 seconds before the ATT moving in and  20 second after ATT moving out, for both LC and BKG data
         This requires using database
        """
        if not self._data:
            print('Light curve are not ready')
            return None
        bkg_data=self._data['BKG_LC']
        ql_data=self._data['QL_LC']
        try:
            rcr=ql_data['rcr']
        except (TypeError, KeyError):
            print('rcr not exists in QL data:')
            return
        

        if not np.any(ql_data['rcr']):
            print('No need for correction')
            return None

        ql_lcs=np.array(list(ql_data['lcs'].values()))
        bkg_lcs=np.array(list(bkg_data['lcs'].values()))
        ql_times=ql_data['time']

        bkg_times=bkg_data['time']

        sel=np.where(ql_data['rcr']==1)
        rcr_in_times=ql_times[sel]

        #print(rcr_in_times)

        tin,tout=rcr_in_times[0],rcr_in_times[-2]
        #print('sel:',time_utils.unix2utc(tin), time_utils.unix2utc(tin-tgap) )
        #print(ql_times-ql_times[0])
        ql_time_sel=np.where( (ql_times>=tin-tgap)  & (ql_times <= tin))
        #here still correct
        #print('ql_time sel:', ql_time_sel)
        ql_time_before_in=ql_times[ql_time_sel]
        #print("QL before in:", [time_utils.unix2utc(x) for x in ql_time_before_in])
        #this is wrong
        ql_before_in=np.squeeze(ql_lcs[:, ql_time_sel])

        
        #calculated the max counts before moved in

        bkg_time_sel=np.where( (bkg_times>=tin)  & (tin+tgap >=bkg_times))

        bkg_after_in=np.squeeze(bkg_lcs[:, bkg_time_sel])

        bkg_time_after_in=bkg_times[bkg_time_sel]
        
        min_index=np.argmin(bkg_after_in[0])
        min_bkg_time=bkg_time_after_in[min_index]

        
        ql_max_cnts=np.max(ql_before_in,axis=1)
        #ql_max_cnts before inserted
        #this works but below is better
        #ql_max_cnts=np.zeros(5)
        #for i in range(5):
        #    ql_max_cnts[i]= extrapolate_counts(ql_time_befor_in, ql_before_in[i], min_bkg_time)
            #estimate  the counts, when bkg counts is the minimal




        bkg_min_cnts=np.min(bkg_after_in,axis=1)

        stat = mdb.get_nearest_lc_stats(ql_time_before_in[0], max_limit=500)

        BKG_baseline=[14,3,3,42,20] #we assume the BKG BKG is stable
        try:
            QL_baseline=stat['median'] #we assume the BKG BKG is stable
        except Exception as e:
            QL_baseline=[271,41,67,799,450] #default background

        slopes=[]
        for i in range(5):

            if bkg_min_cnts[i] <= BKG_baseline[i] + 3* np.sqrt(BKG_baseline[i]) :
                #too noisy, not valid
                slopes.append(0)
            else:
                bkg_counts=bkg_min_cnts[i]-BKG_baseline[i]
                slopes.append((ql_max_cnts[i]-QL_baseline[i])/bkg_counts)


        sel_t_bkg=np.where((bkg_times >= tin) & (tout >=bkg_times))
        bkg_att_in_lcs = []
        for i in range(5):
            bkg_att_in_lcs.append(np.array(bkg_lcs[i][sel_t_bkg]*slopes[i]+QL_baseline[i]) if slopes[i]>0 else bkg_lcs[0][sel_t_bkg]*0)

        bkg_att_in_times=bkg_times[sel_t_bkg]

        return {'time':bkg_att_in_times, 
                'counts': bkg_att_in_lcs, 
                'ql_time_before_in': [time_utils.unix2utc(x) for x in ql_time_before_in],
                'ql_max_cnts_before_in': ql_max_cnts, 
                'bkg_min_cnts_after_in': bkg_min_cnts, 
                'ql_bkg':QL_baseline, 
                'ql_counts_before_in':ql_before_in,
                'BKG_bkg': BKG_baseline,
                'bkg_stat':stat
                }
    





    @classmethod
    def from_database(cls, start_utc, end_utc):
        qlbkg=BKGLigthCurveAnalyzer.from_database(start_utc, end_utc)
        qllc=LightCurveAnalyzer.from_database(start_utc,end_utc)
        return cls(qllc, qlbkg)



