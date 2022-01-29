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
        print('hex:')
        raw = self.packets[self.current_row]['bin']
        text = ''
        for i, h in enumerate(raw):
            text += '{:02X}  '.format(h)
            if i > 0 and (i + 1) % 16 == 0:
                text += '\n'
        print(text)
