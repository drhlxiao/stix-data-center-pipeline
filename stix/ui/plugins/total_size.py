#!/usr/bin/python3
#plugin template
class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row

    def run(self):
        # your code goes here
        print('current row: {}'.format(self.current_row))
        print('Number of packets {}:'.format(len(self.packets)))
        print(len(self.packets))

        counter = {}

        total_size = 0
        for packet in self.packets:
            header = packet['header']
            total_size += header['raw_length']
            apid = header['apid']
            if apid not in counter:
                counter[apid] = 0
            counter[apid] += header['raw_length']
        print('Total size:', total_size)

        print('APID, Length')
        for key, value in counter.items():
            print('{}, {}'.format(key, value))
