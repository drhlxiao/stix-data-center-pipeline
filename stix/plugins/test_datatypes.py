#test stix_datetypes
import os

import sys
sys.path.append(os.path.abspath('../../'))
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt
from stix.core import stix_datatypes as sdt
from stix.spice import stix_datetime


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.iql = 0

    def run(self):
        packet=self.packets[1903]
        pkt= sdt.Packet(packet)
        s=pkt.get_one('NIX00272')
        print(s)


