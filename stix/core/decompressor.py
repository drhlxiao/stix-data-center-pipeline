#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : decompressor.py
# @description:
#               decompression of compressed parameters
from stix.core import logger
from stix.core import config
from stix.core import mongo_db
from stix.core import datatypes
from stix.spice import time_utils as sdt
import numpy as np

logger = logger.get_logger()



MAX_STORED_INTEGER = 1e8
#numbers greater than this value will be converted to float type

DATA_REQUEST_REPORT_SPIDS = [54114, 54115, 54116, 54117, 54143]
QL_REPORT_SPIDS = [54118, 54119, 54121, 54120, 54122]

SKM_GROUPS = {
    'EACC': ("NIXD0007", "NIXD0008", "NIXD0009"),
    'ETRIG': ("NIXD0010", "NIXD0011", "NIXD0012"),
    'LC': ("NIXD0101", "NIXD0102", "NIXD0103"),
    'TriggerSSID30': ("NIXD0104", "NIXD0105", "NIXD0106"),
    'BKG': ("NIXD0108", "NIXD0109", "NIXD0110"),
    'TRIG': ("NIXD0112", "NIXD0113", "NIXD0114"),
    'SPEC': ("NIXD0115", "NIXD0116", "NIXD0117"),
    'VAR': ("NIXD0118", "NIXD0119", "NIXD0120"),
    'CALI': ("NIXD0126", "NIXD0127", "NIXD0128"), 
}
COMPRESSED_PACKET_SPIDS = [
    54112, 54113, 54114, 54115, 54116, 54117, 54118, 54119, 54120, 54121,
    54123, 54124, 54142, 54143, 54110, 54111
]
SCHEMAS = {
    54120: {
        'SKM_Groups':
        ['SPEC', 'TRIG'],  #tell the decompressor to  capture the parameters
        'parameters': {
            'NIX00452': SKM_GROUPS['SPEC'],  #the SKM parameters used to decompress it
            'NIX00453': SKM_GROUPS['SPEC'],
            'NIX00454': SKM_GROUPS['SPEC'],
            'NIX00455': SKM_GROUPS['SPEC'],
            'NIX00456': SKM_GROUPS['SPEC'],
            'NIX00457': SKM_GROUPS['SPEC'],
            'NIX00458': SKM_GROUPS['SPEC'],
            'NIX00459': SKM_GROUPS['SPEC'],
            'NIX00460': SKM_GROUPS['SPEC'],
            'NIX00461': SKM_GROUPS['SPEC'],
            'NIX00462': SKM_GROUPS['SPEC'],
            'NIX00463': SKM_GROUPS['SPEC'],
            'NIX00464': SKM_GROUPS['SPEC'],
            'NIX00465': SKM_GROUPS['SPEC'],
            'NIX00466': SKM_GROUPS['SPEC'],
            'NIX00467': SKM_GROUPS['SPEC'],
            'NIX00468': SKM_GROUPS['SPEC'],
            'NIX00469': SKM_GROUPS['SPEC'],
            'NIX00470': SKM_GROUPS['SPEC'],
            'NIX00471': SKM_GROUPS['SPEC'],
            'NIX00472': SKM_GROUPS['SPEC'],
            'NIX00473': SKM_GROUPS['SPEC'],
            'NIX00474': SKM_GROUPS['SPEC'],
            'NIX00475': SKM_GROUPS['SPEC'],
            'NIX00476': SKM_GROUPS['SPEC'],
            'NIX00477': SKM_GROUPS['SPEC'],
            'NIX00478': SKM_GROUPS['SPEC'],
            'NIX00479': SKM_GROUPS['SPEC'],
            'NIX00480': SKM_GROUPS['SPEC'],
            'NIX00481': SKM_GROUPS['SPEC'],
            'NIX00482': SKM_GROUPS['SPEC'],
            'NIX00483': SKM_GROUPS['SPEC'],
            'NIX00484': SKM_GROUPS['TRIG']
        }
    },
    54124: {
        'SKM_Groups': ['CALI'],
        'parameters': {
            'NIX00158': SKM_GROUPS['CALI']
        }
    },
    54118: {
        'SKM_Groups': ['LC', 'TriggerSSID30'],
        'parameters': {
            'NIX00272': SKM_GROUPS['LC'],
            'NIX00274': SKM_GROUPS['TriggerSSID30']
        }
    },
    54119: {
        'SKM_Groups': ['BKG', 'TRIG'],
        'parameters': {
            'NIX00278': SKM_GROUPS['BKG'],
            'NIX00274': SKM_GROUPS['TRIG']
        }
    },
    54121: {
        'SKM_Groups': ['VAR'],
        'parameters': {
            'NIX00281': SKM_GROUPS['VAR']
        }
    },
    54110: {
        'SKM_Groups': ['EACC', 'ETRIG'],
        'parameters': {
            'NIX00065': SKM_GROUPS['EACC'],
            'NIX00408': SKM_GROUPS['ETRIG'],
            'NIX00409': SKM_GROUPS['ETRIG'],
            'NIX00410': SKM_GROUPS['ETRIG'],
            'NIX00411': SKM_GROUPS['ETRIG'],
            'NIX00412': SKM_GROUPS['ETRIG'],
            'NIX00413': SKM_GROUPS['ETRIG'],
            'NIX00414': SKM_GROUPS['ETRIG'],
            'NIX00415': SKM_GROUPS['ETRIG'],
            'NIX00416': SKM_GROUPS['ETRIG'],
            'NIX00417': SKM_GROUPS['ETRIG'],
            'NIX00418': SKM_GROUPS['ETRIG'],
            'NIX00419': SKM_GROUPS['ETRIG'],
            'NIX00420': SKM_GROUPS['ETRIG'],
            'NIX00421': SKM_GROUPS['ETRIG'],
            'NIX00422': SKM_GROUPS['ETRIG'],
            'NIX00423': SKM_GROUPS['ETRIG']
        }
    },
    54111: {
        'SKM_Groups': ['EACC', 'ETRIG'],
        'parameters': {
            'NIX00260': SKM_GROUPS['EACC'],
            'NIX00242': SKM_GROUPS['ETRIG'],
            'NIX00243': SKM_GROUPS['ETRIG'],
            'NIX00244': SKM_GROUPS['ETRIG'],
            'NIX00245': SKM_GROUPS['ETRIG'],
            'NIX00246': SKM_GROUPS['ETRIG'],
            'NIX00247': SKM_GROUPS['ETRIG'],
            'NIX00248': SKM_GROUPS['ETRIG'],
            'NIX00249': SKM_GROUPS['ETRIG'],
            'NIX00250': SKM_GROUPS['ETRIG'],
            'NIX00251': SKM_GROUPS['ETRIG'],
            'NIX00252': SKM_GROUPS['ETRIG'],
            'NIX00253': SKM_GROUPS['ETRIG'],
            'NIX00254': SKM_GROUPS['ETRIG'],
            'NIX00255': SKM_GROUPS['ETRIG'],
            'NIX00256': SKM_GROUPS['ETRIG'],
            'NIX00257': SKM_GROUPS['ETRIG']
        }
    },
    54112: {
        'SKM_Groups': ['EACC', 'ETRIG'],
        'parameters': {
            'NIX00260': SKM_GROUPS['EACC'],
            'NIX00242': SKM_GROUPS['ETRIG'],
            'NIX00243': SKM_GROUPS['ETRIG'],
            'NIX00244': SKM_GROUPS['ETRIG'],
            'NIX00245': SKM_GROUPS['ETRIG'],
            'NIX00246': SKM_GROUPS['ETRIG'],
            'NIX00247': SKM_GROUPS['ETRIG'],
            'NIX00248': SKM_GROUPS['ETRIG'],
            'NIX00249': SKM_GROUPS['ETRIG'],
            'NIX00250': SKM_GROUPS['ETRIG'],
            'NIX00251': SKM_GROUPS['ETRIG'],
            'NIX00252': SKM_GROUPS['ETRIG'],
            'NIX00253': SKM_GROUPS['ETRIG'],
            'NIX00254': SKM_GROUPS['ETRIG'],
            'NIX00255': SKM_GROUPS['ETRIG'],
            'NIX00256': SKM_GROUPS['ETRIG'],
            'NIX00257': SKM_GROUPS['ETRIG']
        }
    },
    54113: {
        # need to check
        'SKM_Groups': ['EACC', 'ETRIG'],
        'parameters': {
            #'NIX00261': SKM_GROUPS['EACC'],
            'NIX00242': SKM_GROUPS['ETRIG'],
            'NIX00243': SKM_GROUPS['ETRIG'],
            'NIX00244': SKM_GROUPS['ETRIG'],
            'NIX00245': SKM_GROUPS['ETRIG'],
            'NIX00246': SKM_GROUPS['ETRIG'],
            'NIX00247': SKM_GROUPS['ETRIG'],
            'NIX00248': SKM_GROUPS['ETRIG'],
            'NIX00249': SKM_GROUPS['ETRIG'],
            'NIX00250': SKM_GROUPS['ETRIG'],
            'NIX00251': SKM_GROUPS['ETRIG'],
            'NIX00252': SKM_GROUPS['ETRIG'],
            'NIX00253': SKM_GROUPS['ETRIG'],
            'NIX00254': SKM_GROUPS['ETRIG'],
            'NIX00255': SKM_GROUPS['ETRIG'],
            'NIX00256': SKM_GROUPS['ETRIG'],
            'NIX00257': SKM_GROUPS['ETRIG']
        }
    },
    #54142:{},
    54114: {
        # need to check
        'SKM_Groups': ['EACC', 'ETRIG'],
        'parameters': {
            'NIX00065': SKM_GROUPS['EACC'],  #TBC
            'NIX00408': SKM_GROUPS['ETRIG'],
            'NIX00409': SKM_GROUPS['ETRIG'],
            'NIX00410': SKM_GROUPS['ETRIG'],
            'NIX00411': SKM_GROUPS['ETRIG'],
            'NIX00412': SKM_GROUPS['ETRIG'],
            'NIX00413': SKM_GROUPS['ETRIG'],
            'NIX00414': SKM_GROUPS['ETRIG'],
            'NIX00415': SKM_GROUPS['ETRIG'],
            'NIX00416': SKM_GROUPS['ETRIG'],
            'NIX00417': SKM_GROUPS['ETRIG'],
            'NIX00418': SKM_GROUPS['ETRIG'],
            'NIX00419': SKM_GROUPS['ETRIG'],
            'NIX00420': SKM_GROUPS['ETRIG'],
            'NIX00421': SKM_GROUPS['ETRIG'],
            'NIX00422': SKM_GROUPS['ETRIG'],
            'NIX00423': SKM_GROUPS['ETRIG']
        }
    },
    54115: {
        # need to check
        'SKM_Groups': ['EACC', 'ETRIG'],
        'parameters': {
            'NIX00260': SKM_GROUPS['EACC'],  #TBC
            'NIX00242': SKM_GROUPS['ETRIG'],
            'NIX00243': SKM_GROUPS['ETRIG'],
            'NIX00244': SKM_GROUPS['ETRIG'],
            'NIX00245': SKM_GROUPS['ETRIG'],
            'NIX00246': SKM_GROUPS['ETRIG'],
            'NIX00247': SKM_GROUPS['ETRIG'],
            'NIX00248': SKM_GROUPS['ETRIG'],
            'NIX00249': SKM_GROUPS['ETRIG'],
            'NIX00250': SKM_GROUPS['ETRIG'],
            'NIX00251': SKM_GROUPS['ETRIG'],
            'NIX00252': SKM_GROUPS['ETRIG'],
            'NIX00253': SKM_GROUPS['ETRIG'],
            'NIX00254': SKM_GROUPS['ETRIG'],
            'NIX00255': SKM_GROUPS['ETRIG'],
            'NIX00256': SKM_GROUPS['ETRIG'],
            'NIX00257': SKM_GROUPS['ETRIG']
        }
    },
    54116: {
        # need to check
        'SKM_Groups': ['EACC', 'ETRIG'],
        'parameters': {
            'NIX00260': SKM_GROUPS['EACC'],  #TBC
            'NIX00242': SKM_GROUPS['ETRIG'],
            'NIX00243': SKM_GROUPS['ETRIG'],
            'NIX00244': SKM_GROUPS['ETRIG'],
            'NIX00245': SKM_GROUPS['ETRIG'],
            'NIX00246': SKM_GROUPS['ETRIG'],
            'NIX00247': SKM_GROUPS['ETRIG'],
            'NIX00248': SKM_GROUPS['ETRIG'],
            'NIX00249': SKM_GROUPS['ETRIG'],
            'NIX00250': SKM_GROUPS['ETRIG'],
            'NIX00251': SKM_GROUPS['ETRIG'],
            'NIX00252': SKM_GROUPS['ETRIG'],
            'NIX00253': SKM_GROUPS['ETRIG'],
            'NIX00254': SKM_GROUPS['ETRIG'],
            'NIX00255': SKM_GROUPS['ETRIG'],
            'NIX00256': SKM_GROUPS['ETRIG'],
            'NIX00257': SKM_GROUPS['ETRIG']
        }
    },
    54117: {
        # need to check
        'SKM_Groups': ['EACC', 'ETRIG'],
        'parameters': {
            # 'NIX00260': SKM_GROUPS['EACC'],  #TBC
            'NIX00242': SKM_GROUPS['ETRIG'],
            'NIX00243': SKM_GROUPS['ETRIG'],
            'NIX00244': SKM_GROUPS['ETRIG'],
            'NIX00245': SKM_GROUPS['ETRIG'],
            'NIX00246': SKM_GROUPS['ETRIG'],
            'NIX00247': SKM_GROUPS['ETRIG'],
            'NIX00248': SKM_GROUPS['ETRIG'],
            'NIX00249': SKM_GROUPS['ETRIG'],
            'NIX00250': SKM_GROUPS['ETRIG'],
            'NIX00251': SKM_GROUPS['ETRIG'],
            'NIX00252': SKM_GROUPS['ETRIG'],
            'NIX00253': SKM_GROUPS['ETRIG'],
            'NIX00254': SKM_GROUPS['ETRIG'],
            'NIX00255': SKM_GROUPS['ETRIG'],
            'NIX00256': SKM_GROUPS['ETRIG'],
            'NIX00257': SKM_GROUPS['ETRIG']
        }
    },
    54143: {
        # need to check
        'SKM_Groups': ['EACC', 'ETRIG'],
        'parameters': {
            'NIX00268': SKM_GROUPS['EACC'],  #TBC
            'NIX00267': SKM_GROUPS['ETRIG']
        }
    }

    # 54143:{},
}
    

def decompress(x, S, K, M):
    """
    decompress x 
    S, K, M
    """
    if S + K + M > 8 or S not in (0, 1) or K > 7 or M > 7:
        logger.warning('Invalid SKM values: {}{}{}'.format(S, K, M))
        return None
    if K == 0 or M == 0:
        return None

    sign = 1
    if S == 1:  #signed
        MSB = x & (1 << 7)
        if MSB != 0:
            sign = -1
        x = x & ((1 << 7) - 1)

    x0 = 1 << (M + 1)
    if x < x0:
        return x
    mask1 = (1 << M) - 1
    mask2 = (1 << M)
    mantissa1 = x & mask1
    exponent = (x >> M) - 1
    # number of shifted bits
    mantissa2 = mask2 | mantissa1  #add 1 before mantissa
    low = mantissa2 << exponent  #minimal possible value
    high = low | ((1 << exponent) - 1)  #maximal possible value
    mean = (low + high) >> 1  #mean value

    if mean > MAX_STORED_INTEGER:
        return float(mean)

    return sign * mean


class StixDecompressor(object):
    def __init__(self):
        self.is_a_compressed_packet = False
        self.spid = None
        self.SKM_parameters_names = []
        self.SKM_values = dict()
        self.is_a_compressed_packet_parameter_names = []
        self.schema = None
        self.header = {}
        self.header_unix_time=None
        self.unscale_config={'n_trig':1, 'factor':None}
        self.is_trig_scaled_packet=False
        self.mdb= None
        self.packet={}
        self.current_timestamp=None

    def set_mongo_db(self, mdb):
        self.mdb=mdb
    def reset(self):
        self.packet={}
        self.schema = None
        self.is_a_compressed_packet = False
        self.unscale_config={'n_trig':1, 'factor':None}
        self.header = {}
        self.parameters={}
        self.is_trig_scaled_packet=False
        self.current_timestamp=None

    def is_compressed(self):
        return self.is_a_compressed_packet
        
    def add_meta_to_header(self):
        if self.is_trig_scaled_packet and self.header:
            self.header['tsf']=self.unscale_config.get('factor',None)
            self.header['tsfs']=self.unscale_config.get('scaling_factor_source',None)

    def init_decompressor(self):
        
        #self.header=header

        self.is_a_compressed_packet = False
        spid=self.header['SPID']
        self.spid = spid

        if self.spid not in COMPRESSED_PACKET_SPIDS:
            return

        self.is_a_compressed_packet = True

        self.SKM_parameters_names = []
        self.SKM_values = dict()
        self.is_a_compressed_packet_parameter_names = []
        if spid not in SCHEMAS:
            self.is_a_compressed_packet = False
            logger.warning(
                'No schema to decompress packet (SPID {})'.format(
                    spid))
            return
        try:
            self.schema = SCHEMAS[spid]
        except KeyError:
            logger.warning(
                'Failed to decompress packet (SPID {}) due to missing schema'.format(
                    spid))
            self.is_a_compressed_packet = False
            return

        coarse = self.header['coarse_time']
        fine = self.header['fine_time']
        self.header_unix_time = sdt.scet2unix(coarse, fine)


        if self.spid in QL_REPORT_SPIDS:
            self.set_quicklook_scaling_factor(self.header_unix_time)

        SKM_Groups = self.schema['SKM_Groups']
        for grp_name in SKM_Groups:
            self.SKM_parameters_names.extend(SKM_GROUPS[grp_name])
            #list of compressed parameters
    #def get_delta_time(self):
    #    pass


    def capture_config_parameter(self, param_name, raw):
        """
        if parameter is skm, then set skm
        """
        if param_name in self.SKM_parameters_names:
            #is skm, we capture skm value
            self.SKM_values[param_name] = raw
            return True
        elif param_name == 'NIX00405':
            #only applicable to L1, L2, and L3, for L4 data, no integration time
            self.unscale_config['n_int']= int(raw)  #number of integration in units of 0.1 sec
            return True
        elif param_name == 'NIX00407':
            x=int(raw)
            self.unscale_config['n_trig']= sum([(x>>i)&0x1 for i in range(32)])/2
            return True
        elif param_name == 'NIX00404':
            self.current_timestamp=int(raw)
            return True

        elif param_name=='NIX00037':
            uid= int(raw)  #unique id
            self.unscale_config['unique_id'] = uid  #unique id
            factor = None
            source = None
            if self.mdb:
                #load the factor from mdb
                factor=self.mdb.get_sci_trig_scaling_factor(uid)
                source =1 if factor is not None else -1
                try:
                    factor = int(factor) 
                except TypeError:
                    pass
            if factor is None:
                try:
                    factor = int(config.instrument_config['sci_trig_default_scale_factor'])
                    source =2
                except (KeyError, TypeError, ValueError):
                    source = -2


            self.unscale_config['scaling_factor_source']= source
            self.unscale_config['factor']= factor


            logger.info(f"Found scaling factor for {uid}: {factor} {source}")

            if factor==0:
                factor = None

            return True


        return False



    def set_quicklook_scaling_factor(self, unix_time):
        if self.spid in QL_REPORT_SPIDS:
            self.unscale_config['factor']=config.instrument_config['ql_trig_default_scale_factor']
            self.unscale_config['scaling_factor_source']='default'
            #for entry in config.instrument_config['ql_default_scale_factor']:
            #    if unix_time >= entry['time_range'][0] and unix_time<=entry['time_range'][1]:
            #        self.unscale_config['factor']=entry['factors']['ql']

        



    def unscale_triggers(self, scaled_triggers):
        r"""
        Unscale scaled trigger values.

        Trigger values are scaled on board in compression mode 0,0,7 via the following relation

        T_s = T / (factor * n_int * n_trig)

        where `factor` is a configured parameter, `n_int` is the duration in units of 0.1s and
        `n_trig_groups` number of active trigger groups being summed, which depends on the data
        product given by the SSID.

        Parameters
        ----------
        scaled_triggers : int
            Scaled trigger
        """



        try:
            n_group_from_mask = self.unscale_config['n_trig']
            factor = self.unscale_config['factor']
        except Exception:
            n_group_from_mask=1
            print(self.unscale_config)
            raise ValueError(f'No enough information for unscaling triggers!')

        if self.spid!=54143:
            n_int = self.unscale_config['n_int']
            
        if self.spid in [54115, 54114, 54116,54117]:
            #for level0, 1, 2,3, and the group is already 1
            n_group=1
        else:
            n_group=n_group_from_mask
        self.unscale_config['n_trig']=n_group


        self.is_trig_scaled_packet = True
        # Scaled to ints onboard, bins have scaled width of 1, so error is 0.5 times the total factor
        scaling_error = 0.5 * n_group  * n_int * factor if scaled_triggers >0 else 0
        # The FSW essential floors the value so add 0.5 so trigger is the centre of range +/- error
        unscaled_triggers = (scaled_triggers * n_group * n_int * factor) + scaling_error
        #if n_group>=1:


        return unscaled_triggers#, scaling_error**2

    def get_SKM(self, param_name):
        sch=self.schema['parameters']
        if param_name not in sch:
            return None
        try:
            SKM_parameter_names = sch[param_name]
            return (self.SKM_values[SKM_parameter_names[0]],
                    self.SKM_values[SKM_parameter_names[1]],
                    self.SKM_values[SKM_parameter_names[2]])
        except KeyError:
            return None



    def decompress_raw(self, param_name, raw):
        if not self.capture_config_parameter(param_name, raw):
            #if they are not the parameters to be captured, then
            skm = self.get_SKM(param_name)  #compressed raw values
            if not skm:
                return  None
            if skm == (0,0,7):
                try:
                    decom=self.unscale_triggers(raw)
                    print('decompressed:',raw, decom)
                    return decom 
                except Exception as e:
                    return ""
            return decompress(raw, skm[0], skm[1], skm[2])

        return ""


    def decompress_packet(self, packet):
        self.reset()
        if not packet:
            return None
        self.packet=packet
        self.header, self.parameters=self.packet.get('header',{}), self.packet.get('parameters',[])
        self.init_decompressor()

        if not self.is_a_compressed_packet:
            return None
        for param in self.parameters:
            self.walk(param)
        self.add_meta_to_header()
        return True

    def get_packet(self):
        return {'header':self.header, 'parameters':self.parameters}

    def walk(self, node):
        eng = self.decompress_raw(node[0], node[1])
        if eng:
            #don't modify it
            node[2] = eng
        children=node[3]
        #continue
        for child in children:
            self.walk(child)
            





