#plugin example
import pprint
import sys
import pprint
sys.path.append('..')
sys.path.append('.')
from core import stix_packet_analyzer as sta
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
        packet = self.packets[self.current_row]
        analyzer.load(packet)
        #result = analyzer.to_array("NIX00270/NIX00271/*")
        #print(len(result[0][4]))

        dict_ret = analyzer.to_dict()
        pprint.pprint(dict_ret)
