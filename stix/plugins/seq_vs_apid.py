#plugin example
import pprint

from matplotlib import pyplot as plt
from core import stix_packet_analyzer as sta
import sys
analyzer = sta.analyzer()


class Plugin:
    """ don't modify here """

    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        print("Plugin  loaded ...")

    def run(self):
        # your code goes here
        print('current row')
        print(self.current_row)
        #if len(self.packets) > 1:
        #pprint.pprint(self.packets[self.current_row])
        ts = []
        seqs = []
        last_time = -1
        last_seq = -1
        num = analyzer.merge_packets(self.packets, [
            54102,
        ])
        param_values = analyzer.get_merged_parameters()
        print('Spw-link counter:')
        pprint.pprint(param_values['NIXD0079'])
        print('TC20 counter:')
        pprint.pprint(param_values['NIXD0077'])
        print(len(param_values['NIXD0077']))

        num_hk2 = -1
        for packet in self.packets:

            header = packet['header']
            parameters = packet['parameters']

            if header['SPID'] == 54102:
                num_hk2 += 1
                if param_values['NIXD0079'][num_hk2] == 0:
                    continue
                if param_values['NIXD0077'][num_hk2] == 0:
                    continue
            if header['APID'] == 1444:
                if header['time'] > 2771367000:
                    continue
                if header['SPID'] == 54101:
                    continue
                if last_time == -1:
                    last_time = header['time']
                if last_seq == -1:
                    last_seq = header['seq_count'] - 1
                #ts.append(header['time'])
                if header['time'] >= last_time and header[
                        'seq_count'] == last_seq + 1:
                    ts.append(header['time'])
                    seqs.append(header['seq_count'])
                last_time = header['time']
                last_seq = header['seq_count']
        print(ts)

        plt.plot(ts, seqs, 'o')
        plt.show()
