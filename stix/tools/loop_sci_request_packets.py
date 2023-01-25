#this is script is to fixed the issue of Level 1 data, which have multiple pixel mask values in the same packets
# hualin.xiao created on May. 24, 2022

import sys
import copy
import os
import json
from pprint import pprint
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib import pyplot as plt
from stix.core import mongo_db as db
from stix.core import datatypes as sdt
from pathlib import Path
mdb = db.MongoDB()
col_bsd=mdb.get_collection("bsd")
col_pkts=mdb.get_collection("packets")

class L1PacketReader(object):
    '''
        L1,L2 reports analyzer, tasks:
        - reads L1,L2 packets
        - generates synopsis  used by by web pages 
    '''
    def __init__(self):
        self.time = []
        self.pixel_mask = []
        self.first_pixel_mask=None
        self.detector_mask = []
        self.pixel_indexes=None
        self.extract_masks=True
        self.triggers = []
        self.first_mask_children=None
        self.request_id = -1
        self.packet_unix = 0
        self.first_pixel_index =None
        self.pixel_total_counts = [0] * 384



    def get_spectrum_pixel_indexes(self, detector_mask, pixel_mask_array):
        detectors = [0] * 32
        pixel_mask = 0
        for p in pixel_mask_array:
            pixel_mask |= p
        pixels = [0] * 12
        for i in range(32):
            if (detector_mask & (1 << i)) != 0:
                detectors[i] = 1
        for j in range(12):
            if (pixel_mask & (1 << j)) != 0:
                pixels[j] = 1

        pixel_indexes = []
        for i in range(32):
            for j in range(12):
                if detectors[i] == 0 or pixels[j] == 0:
                    continue
                pixel_indexes.append(i * 12 + j)
        return pixel_indexes

    def process_packets(self, packets):
        ipkt=0
        hash_list=[]
        group = {}
        last_time_bin = 0
        current_time = 0
        segs=[0]*4
        num_time_bins=[]
        num_pixel_masks=[]
        last_len=0
        for pkt in packets:
            is_ok = True
            self.request_id = pkt['parameters'][3][1]
            self.packet_unix = pkt['header']['unix_time']
            num_structures = pkt['parameters'][13][1]
            children = pkt['parameters'][13][3]
            for i in range(0, num_structures):
                #print("i:",i)
                offset = i * 22
                self.pixel_mask = [
                        e[1] for e in children[offset + 2][3] if 'NIXG' not in e[0]
                ]

                if not self.pixel_mask:
                    print("Error", pkt['_id'], ' has empty pixel masks')

                this_len=len(self.pixel_mask)
                
                if last_len!=this_len:
                    print("Error", pkt['_id'], ' inconsistent length, this: ', this_len, ' last', last_len )
                
                last_len=this_len


                    









                    



def loop_bsd(bsd_id):
    pkts=mdb.get_packets_of_bsd_request(bsd_id, header_only=False)
    #pkts=mdb.get_packets_by_ids([bsd_id])
    p=L1PacketReader()
    p.process_packets(pkts)
    


if __name__=='__main__':
    if len(sys.argv)!=2:
        print("loop packets <bsd_id>")
    else:
        bsd_id=int(sys.argv[1])
        loop_bsd(bsd_id)

