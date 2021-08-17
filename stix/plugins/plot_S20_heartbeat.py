#plugin example
import sys
from stix.core import stix_packet_analyzer as sta
from matplotlib import pyplot as plt

analyzer = sta.analyzer()

class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row

    
    def run(self):
        # your code goes here
        plt.ion()
        timestamps = []
        seq_counter= []
        heartbeat_values = []
        csv=open('S20.csv','w')
        csv.write('#, timestamp (SCET), seq counter, heartbeat value in S20')
        for i, packet in enumerate(self.packets):
            header = packet['header']
            SCET=float(header['unix_time'])
            seq=header['seq_count']
            timestamps.append(SCET)
            seq_counter.append(seq)
            parameters = packet['parameters']
            heartbeat_value=parameters[9][3][1][1][0]  #NIX00059, 
            SID=parameters[9][3][0][1][0]  #NIX00059, 
            #As can be seen in the GUI, it is in the 9th parameter group. The children is in the 3rd element of the tuple 
            # the 2nd child, the second column is the parameter value. The value is a tuple. 
            heartbeat_values.append(heartbeat_value)
            csv.write('{}, {}, {}, {}\n'.format(i, SCET, seq,heartbeat_value))



        plt.subplot(2,1,1)
        plt.plot(timestamps, heartbeat_values, label='S20 heart heart vs SCET')
        plt.xlabel('Time (s)')
        plt.ylabel('S20 heartbeat value')

        plt.subplot(2,1,2)
        plt.plot(timestamps, seq_counter, label=' sequence counter value vs SCET')
        plt.xlabel('Time (s)')
        plt.ylabel('Seq. counter')
        plt.show()
        
