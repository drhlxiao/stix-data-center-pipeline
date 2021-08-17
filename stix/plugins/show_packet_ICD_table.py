import sys
from core import stix_idb

IDB = stix_idb.stix_idb()


class Plugin:
    """ don't modify here """

    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        print("Plugin loaded ...")

    def run(self):
        # your code goes here
        if len(self.packets) > 1:
            spid = self.packets[self.current_row]['header']['SPID']
            tspd = self.packets[self.current_row]['header']['TPSD']
            print('{:20s}  | {:<80s}  {:10s}  {:10s}'.format(
                'Name', 'Description', 'Position', 'Width'))
            print('-' * 140)
            if tspd > 0:
                ret = IDB.get_variable_packet_structure(spid)
                for row in ret:
                    print('{:20s}  | {:<80s}  {:10d}  {:10d}'.format(
                        row['PCF_NAME'], row['PCF_DESCR'], row['VPD_POS'],
                        row['PCF_WIDTH']))
            else:
                ret = IDB.get_fixed_packet_structure(spid)
                for row in ret:
                    print('{:20s}  | {:<80s}  {:10d}  {:10d}'.format(
                        row['PCF_NAME'], row['PCF_DESCR'], row['PLF_OFFBY'],
                        row['PCF_WIDTH']))
