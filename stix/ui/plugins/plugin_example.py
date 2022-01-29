#plugin example
import pprint


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
        if len(self.packets) > 1:
            pprint.pprint(self.packets[self.current_row])
