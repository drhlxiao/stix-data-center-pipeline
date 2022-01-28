#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : decompressor.py
# @description:
#               decompression of compressed parameters
from stix.core import logger

logger = logger.get_logger()



MAX_STORED_INTEGER = 1e8
#numbers greater than this value will be converted to float type

SKM_GROUPS = {
    'EACC': ("NIXD0007", "NIXD0008", "NIXD0009"),
    'ETRIG': ("NIXD0010", "NIXD0011", "NIXD0012"),
    'LC': ("NIXD0101", "NIXD0102", "NIXD0103"),
    'TriggerSSID30': ("NIXD0104", "NIXD0105", "NIXD0106"),
    'BKG': ("NIXD0108", "NIXD0109", "NIXD0110"),
    'TRIG': ("NIXD0112", "NIXD0113", "NIXD0114"),
    'SPEC': ("NIXD0115", "NIXD0116", "NIXD0117"),
    'VAR': ("NIXD0118", "NIXD0119", "NIXD0120"),
    'CALI': ("NIXD0126", "NIXD0127", "NIXD0128")
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
            'NIX00452':
            SKM_GROUPS['SPEC'],  #the SKM parameters used to decompress it
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
        self.compressed = False
        self.spid = None
        self.SKM_parameters_names = []
        self.SKM_values = dict()
        self.compressed_parameter_names = []
        self.schema = None

    def reset(self):
        self.schema = None
        self.compressed = False

    def is_compressed(self):
        return self.compressed

    def init(self, spid):
        self.compressed = False
        self.spid = spid
        if self.spid not in COMPRESSED_PACKET_SPIDS:
            return
        self.compressed = True
        self.SKM_parameters_names = []
        self.SKM_values = dict()
        self.compressed_parameter_names = []
        if spid not in SCHEMAS:
            self.compressed = False
            logger.warning(
                'A compressed packet (SPID {}) is not decompressed'.format(
                    spid))
            return
        try:
            self.schema = SCHEMAS[spid]
        except KeyError:
            logger.warning(
                'A compressed packet (SPID {}) is not decompressed'.format(
                    spid))
            self.compressed = False
            return

        SKM_Groups = self.schema['SKM_Groups']
        for grp_name in SKM_Groups:
            self.SKM_parameters_names.extend(SKM_GROUPS[grp_name])

    def set_SKM(self, param_name, raw):
        if param_name in self.SKM_parameters_names:
            self.SKM_values[param_name] = raw
            return True
        return False

    def get_SKM(self, param_name):
        if param_name not in self.schema['parameters']:
            return None
        try:
            SKM_parameter_names = self.schema['parameters'][param_name]
            return (self.SKM_values[SKM_parameter_names[0]],
                    self.SKM_values[SKM_parameter_names[1]],
                    self.SKM_values[SKM_parameter_names[2]])
        except KeyError:
            return None

    def decompress_raw(self, param_name, raw):

        if not self.compressed:
            return None
        if not self.set_SKM(param_name, raw):  #they are not  SKM
            skm = self.get_SKM(param_name)  #compressed raw values
            if skm:
                return decompress(raw, skm[0], skm[1], skm[2])
        return None
