import sys
import pprint
from core import packet_analyzer as sta
analyzer = sta.analyzer()

SPID = 54124


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row

    def run(self):
        timestamp = []
        spectra = []
        print('searching for calibration packets')
        num = 0
        for packet in self.packets:

            spid = 0
            try:
                spid = int(packet['header']['SPID'])
            except ValueError:
                #TC
                continue
            if spid != SPID:
                continue

            analyzer.load(packet)
            detector_ids = analyzer.to_array('NIX00159/NIXD0155')[0]
            pixels_ids = analyzer.to_array('NIX00159/NIXD0156')[0]
            spectra = analyzer.to_array('NIX00159/NIX00146/*')[0]
            for i, spec in enumerate(spectra):
                if sum(spec) > 0:
                    num += 1
                    print('Detector %d Pixel %d, counts: %d ' %
                          (detector_ids[i] + 1, pixels_ids[i], sum(spec)))

        print('Total number of non-empty spectra:%d' % num)


def test(filename='../GU/raw/asw164_calibration_test1.dat'):
    from core import stix
    parser = stix.StixTCTMParser()
    parser.parse_file(filename)
    packets = parser.get_decoded_packets()
    p = Plugin(packets)
    p.run()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        test(sys.argv[1])
    else:
        test()
