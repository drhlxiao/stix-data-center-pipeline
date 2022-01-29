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

        total_length = {}
        total_length_SPID={}
        packets_counter={}
        apid_counter={}
        APIDS={}
        DESC={}
        tc_raw_length=0

        total_size = 0
        num_TC = 0
        num_TM = 0
        for packet in self.packets:
            header = packet['header']

            if header['TMTC'] == 'TC':
                if header['name'] != 'ZIX20128':
                    num_TC += 1
                    tc_raw_length+=header['raw_length']
                continue
            if header['SPID'] == 54103:
                #hk4 should not count
                continue
            total_size += header['raw_length']
            num_TM += 1
            apid = header['apid']
            if apid not in total_length:
                total_length[apid] = 0
            total_length[apid] += header['raw_length']
            spid=header['SPID']
            if spid not in total_length_SPID:
            	total_length_SPID[spid]=0
            total_length_SPID[spid]+=header['raw_length']


            if spid not in packets_counter:
            	packets_counter[spid]=0
            packets_counter[spid]+=1

            if apid not in apid_counter:
            	apid_counter[apid]=0
            apid_counter[apid]+=1

            if spid not in APIDS:
                APIDS[spid]=apid
            if spid not in DESC:
                DESC[spid]=header['descr']
        
            
        print('TM Total size:', total_size)
        print('TC total size:', tc_raw_length)
        print('APID, length')
        for key, value in total_length.items():
            print('{},{}, {}'.format(key, apid_counter[key], value))
        print('Number of TM', num_TM)
        print('Number of TC', num_TC)
        print('SPID, length')
        for key, value in total_length_SPID.items():
            print('{}, {},{}, {},{}'.format(key, APIDS[key],  DESC[key],packets_counter[key], value))
