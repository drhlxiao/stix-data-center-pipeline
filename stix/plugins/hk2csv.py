# convert hk2 to csv
# Author Hualin Xiao

import os

import sys
sys.path.append(os.path.abspath('../../'))
from pprint import pprint
import numpy as np
#from matplotlib import pyplot as plt
#from stix.core import stix_datatypes as sdt
#from stix.spice import stix_datetime
import csv

SPID = 54102



class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.iql = 0

    def run(self):
        #packet = sdt.Packet(self.packets[0])
        rows=[]

        #f=open("/home/xiaohl/Desktop/output.csv", "wb")
        #parameter_names=[packet[i]['name'] for i in range(0,93)]
        packet=self.packets[28]
        names=[x[0] for x in packet['parameters'] if 'NIXG' not in x[0] ]
        packet=self.packets[63]
        names.extend([x[0] for x in packet['parameters'] if 'NIXG' not in x[0]])
        print(names)
        print(len(names))
        #f.write(','.join(parameter_names))
        #f.write('\n')
        #for pkt in self.packets:
        #    packet = sdt.Packet(pkt)
        #    eng=[str(packet[i].eng) for i in range(0,93)]
        #    print(','.join(eng))
        #f.write(','.join(parameter_eng))
        #f.write('\n')

        #f.close()


