
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
        hk1=0
        hk2=0

        for packet in self.packets:
            if packet['header']['SPID']==54101:
                hk1+=1
            if packet['header']['SPID']==54102:
                hk2+=1
        print(hk1,hk2)
