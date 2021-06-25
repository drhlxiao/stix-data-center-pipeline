import os
import sys
import array
sys.path.append(os.path.abspath(__file__ + "/../../"))

from core import stix_packet_analyzer as sta
analyzer = sta.analyzer()

from PyQt5 import QtWidgets, QtCore, QtGui

from ROOT import TGraph, TFile, TCanvas, TH1F, gROOT, TBrowser, gSystem, TH2F, gPad


def graph2(x, y, title, xlabel, ylabel):
    n = len(x)
    g = TGraph(n, array.array('d', x), array.array('d', y))
    g.GetXaxis().SetTitle(xlabel)
    g.GetYaxis().SetTitle(ylabel)
    g.SetTitle(title)
    return g


def hist(k, y, title, xlabel, ylabel):
    h2 = TH1F("h%d" % k, "%s; %s; %s" % (title, xlabel, ylabel), 1024, 0, 1024)
    for i, val in enumerate(y):
        h2.SetBinContent(i + 1, val)
    h2.GetXaxis().SetTitle(xlabel)
    h2.GetYaxis().SetTitle(ylabel)
    h2.SetTitle(title)

    return h2


SPID = 54124


class Plugin:
    def __init__(self, packets=[], current_row=0):
        self.packets = packets
        self.current_row = current_row
        self.ical = 0

    def new_folder(self, foldername):
        self.fout.cd()
        self.fout.mkdir(foldername)
        self.fout.cd(foldername)
        self.hcounts = TH2F("hcounts",
                            "Channel counts; Detector; Pixel; Counts ", 32, 0,
                            32, 12, 0, 12)

    def run(self, filename, pdf=None):
        print('Number of packets : {}'.format(len(self.packets)))
        dirname = os.path.dirname(filename)
        if not filename.endswith('.root'):
            print('Invalid filename')
            return

        self.fout = TFile(filename, 'recreate')
        self.fout.cd()

        spectra_container = []
        T0 = 0
        for packet in self.packets:
            try:
                if int(packet['header']['SPID']) != SPID:
                    continue
            except ValueError:
                continue
            header = packet['header']
            seg_flag = header['seg_flag']

            if seg_flag in (1, 3):
                self.new_folder('calibration_{}'.format(self.ical))

            analyzer.load_packet(packet)
            detector_ids = analyzer.to_array('NIX00159/NIXD0155')[0]
            pixels_ids = analyzer.to_array('NIX00159/NIXD0156')[0]
            spectra = analyzer.to_array('NIX00159/NIX00146/*')[0]
            for i, spec in enumerate(spectra):
                if sum(spec) > 0:
                    det = detector_ids[i]
                    pixel = pixels_ids[i]
                    print('Detector %d Pixel %d, counts: %d ' %
                          (det, pixel, sum(spec)))

                    xlabel = 'Energy channel'
                    ylabel = 'Counts'
                    title = ('Detector %d Pixel %d ' % (det + 1, pixel))
                    g = hist(i, spec, title, xlabel, ylabel)
                    g.SetName('Det{}_Pixel{}'.format(det, pixel))
                    g.SetTitle('Det{}_Pixel{}'.format(det, pixel))
                    spectra_container.append((det, pixel, g))
                    self.hcounts.Fill(det + 0.5, pixel + 0.5, sum(spec))

                    #current_idx+=1
            if seg_flag == 2:
                num = len(spectra_container)
                if num > 0:
                    self.hcounts.Write('hcounts')

                canvas = []
                for i in range(0, 12):
                    c = TCanvas()
                    canvas.append(c)
                for spec in spectra_container:
                    det = spec[0]
                    pixel = spec[1]
                    canvas[pixel].cd()
                    spec[2].SetLineColor(det + 1)
                    spec[2].Draw("same")
                    #spec[2].Write()

                for i in range(0, 12):
                    canvas[i].cd()
                    gPad.BuildLegend()
                    canvas[i].Write()
                spectra_container.clear()
                self.ical += 1
        self.fout.Close()

    #print('spectra saved to calibration.root')
    #print('Total number of non-empty spectra:%d' % num)
