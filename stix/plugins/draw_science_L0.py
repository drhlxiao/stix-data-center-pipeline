# -*- encoding: utf-8 -*-
""" A Qt GuI plugin to plot bulk science level 0 data 

In order to use it, you need to select a bulk science level 0 report in the 
GUI and execute this plugin in the plugin manager
"""
import os
import sys
import numpy as np

from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from datetime import datetime
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from stix.core import stix_datatypes as sdt
SPID = 54114

class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QGridLayout(self._main)

        self._canvas_counts = FigureCanvas(Figure(figsize=(9, 5)))
        layout.addWidget(self._canvas_counts,0,0)
        self.addToolBar(QtCore.Qt.BottomToolBarArea,
                        NavigationToolbar(self._canvas_counts, self))

        self._fig = self._canvas_counts.figure
        self.duration=0

        self._ax_counts = self._canvas_counts.figure.add_subplot(221)
        self._ax_spectrum =self._canvas_counts.figure.add_subplot(222) 
        self._ax_trig = self._canvas_counts.figure.add_subplot(223) 
        self._ax_rcr= self._canvas_counts.figure.add_subplot(224) 

    def plot(self, start,  triggers, rcr,dt, hitmap):

        self._ax_trig.clear()
        self._ax_rcr.clear()
        for i in range(0,16):
           self._ax_trig.plot(start,triggers[i], label='Accum. '+str(i))
        self._ax_trig.set_title('Accumulator')
        self._ax_trig.set_xlabel('SCET (s)')
        self._ax_trig.set_ylabel('Counts in {} s'.format(dt[0]))
        self.duration=start[-1]-start[0]
        self._ax_trig.legend()
        self._ax_trig.set_title('Accumulators')

        self._ax_rcr.clear()
        self._ax_rcr.plot(start,rcr)
        self._ax_rcr.set_xlabel('SCET (s)')
        self._ax_rcr.set_ylabel('RCR')
        self.hitmap=hitmap
        ix = range(0,33)
        iy = range(0,12)
        z=[]
        x=[]
        y=[]
        for det in ix:
            for pixel in iy:
                x.append(det)
                y.append(pixel)
                z.append(np.sum(hitmap[det][pixel][:]))
        nbins = [33, 12]
        h = self._ax_counts.hist2d(
            x,
            y,
            nbins,
            np.array([(0, 33), (0, 12)]),
            weights=z,
            cmin=1,
            cmap=plt.cm.jet)
        self._ax_counts.set_xticks(range(0, 33, 1))
        self._ax_counts.set_yticks(range(0, 12, 1))
        self._ax_counts.set_xlabel('Detector')
        self._ax_counts.set_ylabel('Pixel')
        self._ax_counts.set_title('Counts accumulated in {} s'.format(self.duration))
        self._fig.colorbar(h[3],ax=self._ax_counts)
        cid = self._fig.canvas.mpl_connect('button_press_event', self.on_click)
        self._fig.canvas.draw()
        self._fig.tight_layout()


    def on_click(self, event):
        self._ax_spectrum.clear()
        x = int(event.xdata)
        y = int(event.ydata)
        try:
            data = self.hitmap[x][y][:]
            self._ax_spectrum.plot(data)
            self._ax_spectrum.set_xlabel('Energy channel')
            self._ax_spectrum.set_title('Spectrum')
            self._ax_spectrum.set_ylabel('Counts in {} s'.format(self.duration))
            self._fig.canvas.draw()
        except Exception as e:
            print(str(e))


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.app = ApplicationWindow()

    def run(self):
        num_packets=len(self.packets)
        start=[]
        rcr=[]
        dt=[]
        pmask=[]
        dmask=[]
        triggers={ x: [] for x in range(0,16) }
        spectra=[]
        hitmap=np.zeros((33,12,32))
        #detector id 0 - 31 or 1 -32 
        

        while self.current_row < num_packets:
            pkt = self.packets[self.current_row]
            packet = sdt.Packet(pkt)
            self.current_row+=1
            if not packet.isa(SPID):
                print('Not a L0 bulk science packet')
                continue
            print("loading packet #", self.current_row)
            T0=packet[12].raw
            num=packet[13].raw
            children=packet[13].children
            accu_S=packet[5].raw
            accu_K=packet[6].raw
            accu_M=packet[7].raw
            trig_S=packet[9].raw
            trig_K=packet[10].raw
            trig_M=packet[11].raw

            trig_idx = 1 
            if trig_K +trig_M>0:
                trig_idx=2

            counts_idx= 1 #raw value
            if accu_K +accu_M>0:
                counts_idx=2


            for i in range(0, num):
                offset=i*23
                k=0
                start.append(children[offset][1]*0.1+T0)
                rcr.append(children[offset+1][1])
                dt.append(children[offset+2][1]*0.1)
                pmask.append(children[offset+4][1])
                dmask.append(children[offset+5][1])
                spectrum=np.zeros((32,12,32))

                for off in range(offset+6,offset+22):
                    triggers[k].append(children[off][trig_idx])
                    k+=1
                num_samples=children[22][1]
                samples=children[22][3]
                for j in range(0,num_samples):
                    offset_2=j*5
                    pixel=samples[offset_2+1][1]
                    detector=samples[offset_2+2][1]
                    energy_bin=samples[offset_2+3][1]
                    num_bits=samples[offset_2+4][1]
                    continous_counts=samples[offset_2+4][3]
                    counts=1
                    if num_bits==1:
                        counts=continous_counts[0][counts_idx]
                    if num_bits==2:
                        counts=continous_counts[0][counts_idx]<<8
                        counts+=continous_counts[1][counts_idx]
                    spectrum[detector][pixel][energy_bin]=counts
                    hitmap[detector][pixel][energy_bin]+=counts
                spectra.append(spectrum)
            if packet.seg_flag in [2,3]:
                #note that there is a bug in ASW v167.
                #seg_flag is wrong for the last packet
                break

        self.app.show()
        self.app.plot(start,triggers,rcr,dt, hitmap)


                
    

                



                
    

                




