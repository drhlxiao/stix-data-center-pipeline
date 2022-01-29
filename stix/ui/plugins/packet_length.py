
#!/usr/bin/python3
#plugin template
class Plugin:
    def __init__(self,  packets=[], current_row=0):
        self.packets=packets
        self.current_row=current_row
    def run(self):
        # your code goes here
        print('current row: {}'.format(self.current_row))
        print('Number of packets {}:'.format(len(self.packets)))
        print(len(self.packets))

        for packet in self.packets:
            header=packet['header']
            length=header['length']-9+16
