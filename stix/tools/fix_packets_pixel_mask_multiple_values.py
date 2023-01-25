#this is script is to fixed the issue of Level 1 data, which have multiple pixel mask values in the same packets
# hualin.xiao created on May. 24, 2022
# Packet example https://datacenter.stix.i4ds.net/view/packet/id/3385844

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
        for pkt in packets:
            is_ok = True
            #packet = sdt.Packet(pkt)
            self.request_id = pkt['parameters'][3][1]
            self.packet_unix = pkt['header']['unix_time']
            num_structures = pkt['parameters'][13][1]
            print('packet id',pkt['_id'])
            

            pkt['original_parameters']=copy.deepcopy(pkt['parameters'])


            children = pkt['parameters'][13][3]
            for i in range(0, num_structures):
                #print("i:",i)
                fixable=False
                offset = i * 22
                self.pixel_mask = [
                        e[1] for e in children[offset + 2][3] if 'NIXG' not in e[0]
                ]

                self.detector_mask = children[offset + 3][1]
                self.pixel_indexes = self.get_spectrum_pixel_indexes(
                self.detector_mask, self.pixel_mask)
                if self.first_pixel_index is None:

                    self.first_pixel_index = copy.deepcopy(self.pixel_indexes)
                    self.first_pixel_mask=copy.deepcopy(self.pixel_mask)
                    print("first pixel mask:", children[offset+2])
                    self.first_mask_children=copy.deepcopy(children[offset + 2])

                if self.pixel_indexes == self.first_pixel_index:
                    #good packets
                    continue
                


                is_ok=False

                print(pkt['_id'], " is a bad packet, pixel mask:")
                print("Fist block pixel indexes:", self.first_pixel_index)
                print('This mask:', self.pixel_indexes)
                
                self.pixel_indexes = self.first_pixel_index
                children[offset + 2]=self.first_mask_children


                num_samples = children[21 + offset][1]
                samples = children[21 + offset][3]
                counts = []
                first_length=len(self.first_pixel_index)

                if num_samples==0:
                    print("Nsample ==0")


                for j in range(0, num_samples):#interate over time bins
                    k = j * 4
                    num_ele=len(samples[k + 3][3])
                    print(num_ele, "<-last, this len-> ",  first_length)
                    if num_ele<first_length:
                        print(f"ERROR: Packet {pkt['_id']} can not be fixed")
                        continue
                    #print(num_ele, " will be reduced to => ",  first_length)
                    new_children=[]
                    for idx, e in enumerate(samples[k + 3][3]):
                        #iterate over energies
                        pixel=self.pixel_indexes[idx] #get pixel id 
                        if pixel  in self.first_pixel_index: #if the pixel id is in 
                            new_children.append(copy.deepcopy(e))
                            #copy element
                    samples[k + 3][3]=new_children
                    samples[k + 3][1]=first_length

                    if len(new_children)==first_length:
                        fixable=True
                    else:
                        print(f"ERROR: Packet {pkt['_id']} can not be fixed")
                        continue







                #print("new children has ", len(new_children), " elements!")
                #modify the length
            if not is_ok and fixable:
                #print("replace children : ", children[offset+2],'=>')
                #print(self.first_mask_children)
                #now pkt is corrected
                print("Fixing packet:", pkt['_id'])
                pkt['comments']='Bad pixel mask, fixed on May 25, 2022'
                col_pkts.replace_one({'_id':pkt['_id']},pkt)
                



def fix_bsd(bsd_id):
    pkts=mdb.get_packets_of_bsd_request(bsd_id, header_only=False)
    p=L1PacketReader()
    p.process_packets(pkts)
    


if __name__=='__main__':
    if len(sys.argv)!=2:
        print("fix_packet <bsd_id>")
    else:
        bsd_id=int(sys.argv[1])
        fix_bsd(bsd_id)

