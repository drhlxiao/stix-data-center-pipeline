#plugin example
import sys
from stix.core import datatypes as sdt
from matplotlib import pyplot as plt


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row

    def run(self):
        # your code goes here
        plt.ion()
        timestamp = []
        A0_V = []
        A1_V = []
        B0_V = []
        B1_V = []
        for pkt in self.packets:
            packet = sdt.Packet(pkt)
            if not packet.isa(54102):
                continue
            timestamp.append(float(packet['unix_time']))
            names = ['NIX00078', 'NIX00079', 'NIX00080', 'NIX00081']
            A0_V.append(packet[86]['raw'])
            A1_V.append(packet[87]['raw'])
            B0_V.append(packet[88]['raw'])
            B1_V.append(packet[89]['raw'])

        plt.plot(timestamp, A0_V, label='A0_V')
        plt.plot(timestamp, A1_V, label='A1_V')
        plt.plot(timestamp, B0_V, label='B0_V')
        plt.plot(timestamp, B1_V, label='B1_V')
        plt.xlabel('Time (s)')
        plt.ylabel('Raw voltage ')
        plt.legend()
        plt.show()
        plt.savefig('aspect.png')
