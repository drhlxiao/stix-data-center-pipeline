#this is script is to fixed the issue of Level 1 data, which have multiple pixel mask values in the same packets
# hualin.xiao created on May. 24, 2022
# used to fixed the issue reported at https://github.com/i4Ds/STIX-FSW/issues/1030
# this script force the pixel mask to 0xFF



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
from stix.fits import fits_creator as fc
from stix.analysis import sci_packets_analyzer as sp
from pathlib import Path
mdb = db.MongoDB()
col_bsd=mdb.get_collection("bsd")
col_pkts=mdb.get_collection("packets")
good_pixel_mask=[ 
                    "NIX00442", 
                    8, 
                    "", 
                    [ 
                        [ 
                            "NIXG0404", 
                            1, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXD0407", 
                            1, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXG0404", 
                            2, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXD0407", 
                            2, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXG0404", 
                            4, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXD0407", 
                            4, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXG0404", 
                            8, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXD0407", 
                            8, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXG0404", 
                            16, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXD0407", 
                            16, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXG0404", 
                            32, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXD0407", 
                            32, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXG0404", 
                            64, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXD0407", 
                            64, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXG0404", 
                            128, 
                            "", 
                            []
                        ], 
                        [ 
                            "NIXD0407", 
                            128, 
                            "", 
                            []
                        ]
                    ]
                ] 
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

    def process_packets(self, packets, is_recovery=False):
        ipkt=0
        hash_list=[]
        group = {}
        last_time_bin = 0
        current_time = 0
        segs=[0]*4
        for pkt in packets:
            is_ok = True
            #packet = sdt.Packet(pkt)
            if is_recovery:
                print('recovering: ',pkt['_id'])
                try:
                    pkt['parameters']=copy.deepcopy(pkt['original_parameters'])
                except KeyError:
                    pass

                col_pkts.replace_one({'_id':pkt['_id']},pkt)
                continue
            print('fixing: ',pkt['_id'])

            self.request_id = pkt['parameters'][3][1]
            self.packet_unix = pkt['header']['unix_time']
            num_structures = pkt['parameters'][13][1]
            

            pkt['original_parameters']=copy.deepcopy(pkt['parameters'])
            fixable=True


            children = pkt['parameters'][13][3]
            for i in range(0, num_structures):
                #print("i:",i)
                offset = i * 22
                self.pixel_mask = [
                        e[1] for e in children[offset + 2][3] if 'NIXG' not in e[0]
                ]

                self.detector_mask = children[offset + 3][1]
                self.pixel_indexes = self.get_spectrum_pixel_indexes(
                self.detector_mask, self.pixel_mask)
                if self.first_pixel_index is None:

                    self.first_pixel_index =  self.get_spectrum_pixel_indexes(
                                self.detector_mask, [1, 2, 4, 8, 16, 32, 64, 128])
                    print("first pixel mask:", children[offset+2][0:8] )
                    #self.first_mask_children=copy.deepcopy(children[offset + 2][0:8])
                    #the first 8 element , need to change if the mask is not 0xFF or 0xFFF
                    first_length=len(self.first_pixel_index)

                if self.pixel_indexes == self.first_pixel_index:
                    #good packets
                    continue
                


                is_ok=False

                print(pkt['_id'], " is a bad packet, pixel mask:")
                #print("Fist block pixel indexes:", self.first_pixel_index)
                #print('This mask:', self.pixel_indexes)
                
                pkt['parameters'][13][3][offset + 2]=good_pixel_mask
                #correct the pixel mask



                num_samples = children[21 + offset][1]
                samples = children[21 + offset][3]
                counts = []

                if num_samples==0:
                    print("No sample ==0")
                    continue


                for j in range(0, num_samples):#interate over time bins
                    k = j * 4
                    num_ele=len(samples[k + 3][3])
                    print(num_ele, "<-last, this len-> ",  first_length)
                    if num_ele < first_length:
                        print(f"ERROR: Packet {pkt['_id']} can not be fixed: {num_ele} < {first_length}")
                        continue
                    #print(num_ele, " will be reduced to => ",  first_length)
                    new_children=[]
                    for idx, e in enumerate(samples[k + 3][3]):
                        #iterate over energies
                        pixel=self.pixel_indexes[idx] #get pixel id 

                        if pixel  in self.first_pixel_index: #if the pixel id is in 
                            new_children.append(copy.deepcopy(e))
                            #copy element
                    pkt['parameters'][13][3][21+offset][3][k + 3][3]=copy.deepcopy(new_children)
                    pkt['parameters'][13][3][21+offset][3][k + 3][1]= first_length

                    if len(new_children)==first_length:
                        #fixable=True
                        print("Fixable")
                    else:
                        fixable=False
                        print(f"ERROR: Packet {pkt['_id']} can not be fixed, because the length inconstent: this={len(new_children)}, first={first_length}")
                        print("first pixel mask:", self.first_pixel_index)
                        print("this pixel mask:", self.pixel_indexes)
                        print("sample len:", len(samples[k+3][3]))
                        continue

                self.pixel_indexes = self.first_pixel_index






            #print("new children has ", len(new_children), " elements!")
            #modify the length
            if not is_ok and fixable:
                #print("replace children : ", children[offset+2],'=>')
                #print(self.first_mask_children)
                #now pkt is corrected
                print("Fixing packet:", pkt['_id'])
                pkt['comments']=f'Bad pixel mask, fixed on {datetime.now()}'
                col_pkts.replace_one({'_id':pkt['_id']},pkt)
            else:
                print("not fixable!", is_ok, fixable)
                    



def fix_bsd(bsd_id, is_recovery=False):
    print("Correcting packets for bsd ", bsd_id)
    pkts=mdb.get_packets_of_bsd_request(bsd_id, header_only=False)
    p=L1PacketReader()
    p.process_packets(pkts, is_recovery)
    


if __name__=='__main__':
    if len(sys.argv)==1:
        print("fix_packet <bsd_id> <bsd_id 2 > ...")
    else:
        bsd_ids=[  int(x) for x  in sys.argv[1:]]
        for _id in bsd_ids:

            print("recovery the backup : ", _id)
            fix_bsd(_id, is_recovery=True)
            print("fixing id: ", _id)
            fix_bsd(_id, is_recovery=False)
            print('creating fits ',_id)
            path='/data/fits/'
            fc.create_fits_for_bulk_science(_id, _id, path, overwrite=True, version=1)
            print('creating json: ',_id)
            sp.process_science_request(_id)

