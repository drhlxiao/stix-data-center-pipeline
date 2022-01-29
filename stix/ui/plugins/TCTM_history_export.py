#plugin example
import pprint

ISTART = 0
ISTOP = 0


class Plugin:
    """ don't modify here """

    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        print("Plugin  loaded ...")

    def run(self, istart=ISTART, istop=ISTOP):
        last_time = 0
        num_TC = 0
        total_size = 0
        num_TM = 0
        TC_size = 0
        counter = {}
        packets = self.packets
        if istop > istart:
            packets = self.packets[istart:istop]

        for packet in packets:
            header = packet['header']
            if header['service_type'] == 20 and header[
                    'service_subtype'] == 128:
                continue
            leng = header['raw_length']
            spid = header['SPID']
            extra = ''
            extra2 = ''
            if header['TMTC'] == 'TC':
                if header['name'] == 'ZIX20128':
                    continue
                num_TC += 1
                TC_size += header['raw_length']
                leng = 0
                extra2 = ';'
            else:
                if header['SPID'] in [54103, 54101]:
                    continue
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

            #print('{}; {} {}({},{}) - {}; {} {}; {}'.format(header['UTC'], extra, header['TMTC'], int(header['service_type']),
            #                                      int(header['service_subtype']),
            #                                      header['descr'], extra2, leng, header['apid'],))

        print('Number of TC,', num_TC)
        print('TC size,', TC_size)
        print('TM Total size,', total_size)
        print('Number of TM,', num_TM)
        print('APID, length')
        for key, value in counter.items():
            print('{}, {}'.format(key, value))
