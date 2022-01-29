# -*- encoding: utf-8 -*-
""" A  Qt parser GUI plugin used to plot ASIC on-demand readouts, threshold readouts and 
    dark current readouts. 

In order to use it, you need to select a background report in the 
GUI and execute this plugin in the plugin manager
"""

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
from stix.core import datatypes as sdt

SPIDs = [54132, 54133, 54134, 54130]
#On demand readout or baseline

DESCR = {
    54132: 'On-demand readout',
    54133: 'Baseline readouts',
    54134: 'Dark current readouts',
    54130: 'threshold scan'
}
PIXEL_TO_ASIC= {
    0: 26,
    1: 15,
    2: 8,
    3: 1,
    4: 29,
    5: 18,
    6: 5,
    7: 0,
    8: 30,
    9: 21,
    10: 11,
    11: 3
}
#PIXEL: ASICID
BIG_PIXELS= [26,15,8,1,29,18,5,0]
SMALL_PIXELS= [30,21, 11,3]


def draw_thresholds(description, UTC, detectors, channels, ADC_mean):
    x = np.array(detectors)
    y = np.array(channels)
    chan = x * 32 + y
    z1 = np.array(ADC_mean)

    big_pixel_z=[]
    small_pixel_z=[]
    for ch, adc in zip(channels, ADC_mean):
        if ch in SMALL_PIXELS:
            small_pixel_z.append(adc)
        elif ch in BIG_PIXELS:
            big_pixel_z.append(adc)



    fig = plt.figure()
    fig.suptitle('{} at {}'.format(description, UTC))
    ax = plt.subplot(311)
    nbins = [32, 12]
    h = plt.hist2d(
        x,
        y,
        nbins,
        np.array([(1, 33), (0, 12)]),
        weights=z1,
        cmin=1,
        cmap=plt.cm.jet)
    ax.set_xticks(range(1, 33, 2))
    ax.set_yticks(range(0, 12, 1))
    plt.title('ADC')
    plt.xlabel('Detector #')
    plt.ylabel('Channel ')
    plt.colorbar(h[3], ax=ax)

    ax = plt.subplot(312)
    plt.plot(chan, z1)
    plt.xlabel('channel #')
    plt.ylabel('ADC mean')

    plt.subplot(313)
    plt.hist(z1, bins=100, label='all')
    plt.hist(small_pixel_z,label='small pixels')
    plt.hist(big_pixel_z, label='big pixels')
    plt.xlabel('ADC channel')
    plt.ylabel('Nb. of channels')
    plt.legend()
    fig.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()


def draw_ADC_mean_dev(description, UTC, detectors, channels, ADC_mean,
                      ADC_stdev):
    x = np.array(detectors)
    y = np.array(channels)
    chan = x * 32 + y
    z1 = np.array(ADC_mean)
    z2 = np.array(ADC_stdev)
    big_pixel_z=[]
    small_pixel_z=[]
    for ch, adc in zip(channels, ADC_mean):
        if ch in SMALL_PIXELS:
            small_pixel_z.append(adc)
        elif ch in BIG_PIXELS:
            big_pixel_z.append(adc)


    fig = plt.figure()
    fig.suptitle('{} at {}'.format(description, UTC))
    ax = plt.subplot(211)
    nbins = [32, 12]
    h = plt.hist2d(
        x,
        y,
        nbins,
        np.array([(1, 33), (0, 12)]),
        weights=z1,
        cmin=1,
        cmap=plt.cm.jet)
    ax.set_xticks(range(1, 33, 2))
    ax.set_yticks(range(0, 12, 1))
    plt.title('ADC')
    plt.xlabel('Detector #')
    plt.ylabel('Channel ')
    plt.colorbar(h[3], ax=ax)

    ax = plt.subplot(212)
    h2 = plt.hist2d(
        x,
        y,
        nbins,
        np.array([(1, 33), (0, 12)]),
        weights=z2,
        cmin=1,
        cmap=plt.cm.jet)
    ax.set_xticks(range(1, 33, 2))
    ax.set_yticks(range(0, 12, 1))
    plt.colorbar(h2[3], ax=ax)
    plt.title('ADC std dev.')
    plt.xlabel('Detector #')
    plt.ylabel('Channel ')

    fig.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()

    fig = plt.figure()
    fig.suptitle('{} at {}'.format(description, UTC))
    plt.subplot(221)
    plt.plot(chan, z1)
    plt.xlabel('channel #')
    plt.ylabel('ADC mean')
    plt.subplot(222)
    plt.plot(chan, z2)
    plt.xlabel('channel #')
    plt.ylabel('ADC std. dev.')

    plt.subplot(223)
    plt.hist([z1, small_pixel_z,big_pixel_z], bins=100, histtype='step', label=['all','small','big'])
    #plt.hist(small_pixel_z, histtype='step',label='small')
    #plt.hist(big_pixel_z, histtype='step',label='big')
    plt.legend()
    plt.xlabel('ADC channel')
    plt.ylabel('Nb. of channels')
    plt.subplot(224)
    plt.hist(z2, bins=100, histtype='step')
    plt.xlabel('ADC channel dev.')
    plt.ylabel('Nb. of channels')
    fig.tight_layout()
    plt.subplots_adjust(top=0.85)

    plt.show()


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.iql = 0

    def run(self):
        isub = 0
        fig = None
        detectors = []
        channels = []
        ADC_mean = []
        ADC_stdev = []
        is_thr = False

        num_read = 0

        length = len(self.packets)

        while self.current_row < length:
            pkt = self.packets[self.current_row]
            packet = sdt.Packet(pkt)
            if packet.SPID not in SPIDs and num_read == 0:
                print('Not ADC readout report')
                return

            num_read += 1

            if packet.SPID not in SPIDs:
                print('Not ADC readout report')
                continue
            if packet.SPID == 54130:
                is_thr = True

            description = DESCR[packet.SPID]
            UTC = packet['header']['UTC']
            print('Analyzing packet # {}'.format(self.current_row))
            det = packet.get('NIX00104/NIX00100')[0]
            chan = packet.get('NIX00104/NIX00105')[0]
            adc = []
            std = []
            selector = 'NIX00104/NIX00106'
            if is_thr:
                selector = 'NIX00104/NIX00108'
            else:
                std = packet.get('NIX00104/NIX00107')[0]
            adc = packet.get(selector)[0]
            detectors.extend(det)
            channels.extend(chan)
            ADC_mean.extend(adc)
            ADC_stdev.extend(std)
            if packet.seg_flag in [2, 3]:
                break
            self.current_row += 1
        if is_thr:
            draw_thresholds(description, UTC, detectors, channels, ADC_mean)
        else:
            draw_ADC_mean_dev(description, UTC, detectors, channels, ADC_mean,
                              ADC_stdev)
