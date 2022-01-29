# -*- encoding: utf-8 -*-
""" A Qt GuI plugin to plot calibration spectra 

In order to use it, you need to select a calibration report in the 
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
from stix.core import datatypes as sdt
SPID = 54124


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)
        self.h2counter = None
        self.spectra = None
        self._canvas = FigureCanvas(Figure(figsize=(9, 5)))
        layout.addWidget(self._canvas)
        self.addToolBar(QtCore.Qt.BottomToolBarArea,
                        NavigationToolbar(self._canvas, self))

        self._spectrum_canvas = FigureCanvas(Figure(figsize=(9, 5)))
        layout.addWidget(self._spectrum_canvas)
        self.addToolBar(NavigationToolbar(self._spectrum_canvas, self))
        self._fig = self._canvas.figure
        self._ax = self._canvas.figure.subplots()

        self._ax_spectrum = self._spectrum_canvas.figure.subplots()



    def plot(self, h2counter, spectra):
        self.h2counter = h2counter
        self.spectra = spectra

        x = np.array([e[0] for e in self.h2counter])
        y = np.array([e[1] for e in self.h2counter])
        z = np.array([e[2] for e in self.h2counter])
        nbins = [32, 12]
        h = self._ax.hist2d(
            x,
            y,
            nbins,
            np.array([(0, 32), (0, 12)]),
            weights=z,
            cmin=1,
            cmap=plt.cm.jet)
        self._ax.set_xticks(range(0, 32, 2))
        self._ax.set_yticks(range(0, 12, 1))

        self._ax.set_xlabel('Detector')
        self._ax.set_ylabel('Pixel')
        self._fig.colorbar(h[3])
        cid = self._fig.canvas.mpl_connect('button_press_event', self.on_click)
        self._ax.figure.canvas.draw()

    def on_click(self, event):
        self._ax_spectrum.clear()
        x = int(event.xdata)
        y = int(event.ydata)
        try:
            data = self.spectra[(x, y)]
            self._ax_spectrum.plot(data)
            self._ax_spectrum.figure.canvas.draw()
        except Exception as e:
            print(str(e))


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.ical = 0
        self.h2counter = []
        self.spectra_container = {}
        self.app = ApplicationWindow()

    def run(self):
        print('Number of packets : {}'.format(len(self.packets)))
        num_packets = len(self.packets)
        num_read = 0
        while self.current_row < num_packets:
            pkt = self.packets[self.current_row]
            packet = sdt.Packet(pkt)
            self.current_row += 1

            if not packet.isa(SPID):
                continue
            seg_flag = packet['seg_flag']

            detector_ids = packet.get('NIX00159/NIXD0155')[0]
            pixels_ids = packet.get('NIX00159/NIXD0156')[0]
            spectra = packet.get('NIX00159/NIX00146/*')[0]

            live_time = packet[4].raw
            quiet_time = packet[3].raw
            num_read += 1

            for i, spec in enumerate(spectra):
                if sum(spec) > 0:
                    det = detector_ids[i]
                    pixel = pixels_ids[i]
                    self.spectra_container[(det, pixel)] = spec
                    self.h2counter.append((det, pixel, sum(spec)))

            if seg_flag in [2, 3]:
                break
        num = len(self.spectra_container)
        if num == 0:
            print('spectra empty')
            return
        print('Number of packets read:', num_read)
        self.app.show()
        self.app.plot(self.h2counter, self.spectra_container)
