import os
import xml.etree.ElementTree as Et
from collections import defaultdict

from pathlib import Path
from time import perf_counter
from binascii import unhexlify

import numpy as np
from astropy.io import fits
from astropy.table import vstack, unique

from astropy.table.table import Table, QTable
from astropy.time.core import Time
from bitstring import ConstBitStream, ReadError

from stix.io.fits.processors import FitsL0Processor
from stix.core.stix_logger import get_logger
logger = get_logger(__name__)

tm_header = {
    'version': 'uint:3',
    'type': 'uint:1',
    'dflag': 'uint:1',
    'apid': 'uint:7',
    'apid_cat': 'uint:4',
    'seq_group': 'uint:2',
    'seq_count': 'uint:14',
    'data_len': 'uint:16',
    'spare1': 'uint:1',
    'pus_ver': 'uint:3',
    'spare2': 'uint:4',
    'service_type': 'uint:8',
    'service_subtype': 'uint:8',
    'dest_id': 'uint:8',
    'scet_coarse': 'uint:32',
    'scet_fine': 'uint:16'
}

SECONDS_IN_DAY = 24*60*60

class level0:
    def __init__(self, control, data):
        self.level = 'L0'
        self.type = f"{control['service_type'][0]}-{control['service_subtype'][0]}"
        self.control = control
        self.data = data

        synced_times = (self.control['scet_coarse'] >> 31) != 1
        self.control['scet_coarse'] = np.where(synced_times, self.control['scet_coarse'],
                                               self.control['scet_coarse'] ^ 2**31)
        self.control['time_sync'] = synced_times
        scet = self.control['scet_coarse'] + self.control['scet_coarse']/2**16
        min_time_ind = scet.argmin()
        max_time_ind = scet.argmax()
        self.obt_beg = f"{self.control['scet_coarse'][min_time_ind]}:" \
                       f"{self.control['scet_fine'][min_time_ind]:05d}"

        self.obt_end = f"{self.control['scet_coarse'][max_time_ind]}:" \
                       f"{self.control['scet_fine'][max_time_ind]:05d}"

    def __repr__(self):
        return f"<Level0 {self.level}, {self.type}, " \
               f"{self.obt_beg}-{self.obt_end}>"

    def __add__(self, other):
        if not isinstance(other, type(self)):
            raise TypeError(f'Products must of same type not {type(self)} and {type(other)}')

        other.control['index'] = other.control['index'] + self.control['index'].max() + 1
        other.data['control_index'] = other.data['control_index'] + self.control['index'].max() + 1

        control = vstack((self.control, other.control))
        control = unique(control, ['scet_coarse', 'scet_fine', 'seq_count'])

        orig_indices = control['index']
        new_index = range(len(control))
        control['index'] = new_index

        data = vstack((self.data, other.data))
        data = data[np.nonzero(orig_indices[:, None] == data['control_index'].data)[1]]
        data['control_index'] = range(len(data))

        if np.abs([((len(data['data'][i]) / 2) - (control['data_len'][i] + 7))
                   for i in range(len(data))]).sum() > 0:
            logger.error('Expected and actual data length do not match')

        return type(self)(control, data)

    @classmethod
    def from_fits(cls, fitspath):
        control = QTable.read(fitspath, hdu='CONTROL')
        data = QTable.read(fitspath, hdu='DATA')
        return cls(control=control, data=data)

    def to_days(self):
        scet = self.control['scet_coarse'] + (self.control['scet_coarse'] / (2**16-1))
        scet = scet[(self.control['scet_coarse'] >> 31) != 1]
        days, frac = np.divmod(scet, SECONDS_IN_DAY)
        days = np.unique(days) * SECONDS_IN_DAY
        for day_start in days:
            day_end = day_start + SECONDS_IN_DAY

            inds = np.argwhere((scet >= day_start) & (scet < day_end))
            i = np.arange(inds.min(), inds.max()+1)
            control = self.control[i]
            control_indices = control['index']
            min_index = control_indices.min()
            control['index'] = control['index'] - min_index
            data = self.data[i]
            data['control_index'] = data['control_index'] - min_index

            yield type(self)(control=control, data=data)


class SciLevel0(level0):
    def __init__(self, control, data):
        super().__init__(control, data)
        self.ssid = self.control['ssid'][0]

    def to_requests(self):
        yield self

def process_tmtc_file(tmfile, basedir):
    fits_processor = FitsL0Processor(basedir)
    tree = Et.parse(tmfile)
    root = tree.getroot()
    packet_data = defaultdict(list)
    for i, node in enumerate(root.iter('Packet')):
        packet_binary = unhexlify(node.text)
        # Not sure why guess and extra moc header
        header, packet_hex = process_tm_packet(packet_binary[76:])
        # if header.get('ssid', -1) == 21:
        key = f"{header['service_type']}-{header['service_subtype']}"
        packet_data[key].append((header, packet_hex))
    for product, packet in packet_data.items():

        #header and hex for the same type
            # packet_data = sorted(packet_data,
            #        key=lambda x: (x[0]['scet_coarse'] + x[0]['scet_fine'] / 2 ** 16, x[0]['seq_count']))

        #difference services
        header, hex_data = zip(*packet)
        control = Table(header)
        control['index'] = np.arange(len(control))
        data = Table()
        data['control_index'] = control['index']
        # data['data'] = binary_data
        data['data'] = hex_data
        if 'ssid' in control.colnames:
            ssids = np.unique(control['ssid'])
            for ssid in ssids:
                index = np.nonzero(control['ssid'] == ssid)
                if len(index[0]) > 0:
                    cur_control = control[index]
                    cur_data = data[index]
                    cur_control['index'] = np.arange(len(cur_control))
                    cur_data['control_index'] = np.arange(len(cur_control))

                    prod = SciLevel0(control=cur_control, data=cur_data)
                    fits_processor.write_fits(prod)
        else:
            prod = level0(control=control, data=data)
            fits_processor.write_fits(prod)

    # return packet_data


def process_tm_packet(binary_packet):
    bits = ConstBitStream(binary_packet)
    header = {k: v for k, v in
              zip(tm_header.keys(), bits.readlist(tm_header.values()))}
    if header['service_type'] in [3, 21] and header['service_subtype'] in [25, 6]:
        try:
            ssid = bits.peek('uint:8')
            header['ssid'] = ssid
        except ReadError as e:
            header['ssid'] = -1

    if header['data_len']+7 != len(bits.hex)/2:
        raise ValueError('Expected and actual packet lentght do not match')

    return header, bits.hex


if __name__ == '__main__':
    tstart = perf_counter()

    # Real data
    raw_tmtc = Path('/Users/shane/Projects/STIX/dataview/data/real')
    tmtc_files = sorted(list(raw_tmtc.glob('*PktTmRaw*.xml')), key=os.path.getctime)
    # tmtc_file = ['/Users/shane/Projects/STIX/dataview/data/real/IX4.DAY2.BatchRequest.PktTmRaw.SOL.0.2020.140.11.55.36.842.yzKB@2020.140.11.55.37.803.1.xml']
    bd = Path('/Users/shane/Projects/STIX/dataview/data/test2')
    for tmtc_file in tmtc_files:
        process_tmtc_file(tmtc_file, basedir=bd)

    tend = perf_counter()
    logger.info('Time taken %f', tend - tstart)
