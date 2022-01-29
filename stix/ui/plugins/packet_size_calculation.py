#plugin example
import pprint


class Plugin:
    """ don't modify here """
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        print("Plugin  loaded ...")
    def run(self, istart=2500, istop=2574):
        print('current row')
        print(self.current_row)

        last_time = 0
        num_TC = 0
        total_size = 0
        num_TM = 0
        counter = {}
        for packet in self.packets[istart:istop]:
            header = packet['header']
            if header['service_type'] == 20 and header[
                    'service_subtype'] == 128:
                continue
            leng = header['raw_length']
            spid = header['SPID']
            extra = ''
            extra2 = ''
            if header['TMTC'] == 'TC':
                num_TC += 1
                leng = 0
                extra2 = ';'
            else:
                #if int(header['SPID']) not in :
                #    continue
                total_size += header['raw_length']
                num_TM += 1
                apid = header['apid']
                if apid not in counter:
                    counter[apid] = 0
                counter[apid] += header['raw_length']
                extra = ';'
            if spid == 54102:
                if header['SCET'] < last_time + 14:
                    continue
                last_time = header['SCET']
        print('TM Total size,', total_size)
        print('APID, length')
        for key, value in counter.items():
            print('{}, {}'.format(key, value))
        print('Number of TM,', num_TM)
        print('Number of TC,', num_TC)
