#plot calibration spectra, whose counts are still compressed

import os

import sys
from stix.core import packet_analyzer as sta
analyzer = sta.analyzer()

from PyQt5 import QtWidgets, QtCore, QtGui

from ROOT import TGraph, TFile, TCanvas, TH1F, gROOT, TBrowser, gSystem


def graph2(x, y, title, xlabel, ylabel):
    n = len(x)
    g = TGraph(n, array('d', x), array('d', y))
    g.GetXaxis().SetTitle(xlabel)
    g.GetYaxis().SetTitle(ylabel)
    g.SetTitle(title)
    return g


def hist(k, y, title, xlabel, ylabel):
    h2 = TH1F("h%d" % k, "%s; %s; %s" % (title, xlabel, ylabel), 1024, 0, 1024)
    for i, val in enumerate(y):
        h2.SetBinContent(i + 1, val)
        #for j in range(val):
        #    h2.Fill(i)
    h2.GetXaxis().SetTitle(xlabel)
    h2.GetYaxis().SetTitle(ylabel)
    h2.SetTitle(title)
    #h2.SetEntries(sum)

    return h2


SPID = 54124


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row

    def run(self):
        # your code goes here
        print('Number of packets : {}'.format(len(self.packets)))
        filename = str(
            QtWidgets.QFileDialog.getSaveFileName(None, "Save to file", '',
                                                  "ROOT(*.root)")[0])
        dirname = os.path.dirname(filename)
        if not filename.endswith('.root'):
            print('Invalid filename')
            return

        fout = TFile(filename, 'recreate')
        fout.cd()
        hcounts = TH1F("hcounts", "Channel counts; Pixel #; Counts", 12 * 32,
                       0, 12 * 32)
        num = 0
        for packet in self.packets:
            try:
                if int(packet['header']['SPID']) != SPID:
                    continue
            except ValueError:
                continue
            analyzer.load(packet)
            detector_ids = analyzer.to_array(
                'NIX00159/NIXD0155', dtype='eng')[0]
            pixels_ids = analyzer.to_array('NIX00159/NIXD0156', dtype='eng')[0]
            spectra = analyzer.to_array('NIX00159/NIX00146/*', dtype='eng')[0]
            for i, spec in enumerate(spectra):
                if sum(spec) > 0:
                    num += 1
                    det = detector_ids[i]
                    pixel = pixels_ids[i]
                    print('Detector %d Pixel %d, counts: %d ' % (det, pixel,
                                                                 sum(spec)))

                    xlabel = 'Energy channel'
                    ylabel = 'Counts'
                    title = ('Detector %d Pixel %d ' % (det + 1, pixel))
                    g = hist(i, spec, title, xlabel, ylabel)
                    #cc.cd(current_idx+1)
                    #g.Draw("hist")
                    g.Write()
                    hcounts.Fill(12 * det + pixel, sum(spec))
                    #current_idx+=1
        if num > 0:
            hcounts.Write('hcounts')
        print('spectra saved to calibration.root')
        fout.Close()
        gROOT.ProcessLine('new TFile("{}")'.format(filename))
        gROOT.ProcessLine('new TBrowser()')
        print('Total number of non-empty spectra:%d' % num)
