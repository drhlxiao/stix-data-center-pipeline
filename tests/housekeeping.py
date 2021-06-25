import os
import sys
import math

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime

#sys.path.append(os.path.abspath(__file__ + "/../../"))

from stix.core import stix_datatypes as sdt
from stix.utils import stix_desc
from stix.core import stix_idb
STIX_IDB = stix_idb.stix_idb()

SPIDs = [54101, 54102]

groups = [{
    'Operation mode': ['NIXD0023']
}, {
    'ASW': ['NIXD0001', 'NIXD0021', 'NIXD0005', 'NIXD0003', 'NIXD0166']
}, {
    'CPU': [
        'NIXD0002',
        'NIXD0167',
        'NIXD0168',
        'NIXD0169',
    ]
}, {
    'Spw-links': ['NIXD0030', 'NIXD0031', 'NIXD0032', 'NIXD0001']
}, {
    'Spw-links': ['NIXD0078', 'NIXD0079', 'NIXD0080', 'NIXD0081']
}, {
    'Telecommands': ['NIXD0077']
}, {
    'IDPU': ['NIXD0004']
}, {
    'FDIR': ['NIX00085']
}, {
    'DPU voltages':
    ['NIXD0027', 'NIXD0028', 'NIXD0029', 'NIXD0030', 'NIXD0031']
}, {
    'DPU voltages': ['NIXD0032', 'NIXD0035', 'NIXD0036', 'NIXD0037']
}, {
    'DPU temperature sensors': ['NIXD0025', 'NIXD0026']
}, {
    'PSU': [
        'NIXD0024',
    ]
}, {
    'LV': ['NIXD0092']
}, {
    'Aspect system power/voltages':
    ['NIXD0086', 'NIXD0087', 'NIXD0038', 'NIXD0039', 'NIXD0048', 'NIXD0049']
}, {
    'Aspect system temperature sensors':
    ['NIXD0040', 'NIXD0041', 'NIXD0042', 'NIXD0043']
}, {
    'Aspect system temperature sensors group 2':
    ['NIXD0044', 'NIXD0045', 'NIXD0046', 'NIXD0047']
}, {
    'Aspect readouts': ['NIX00078', 'NIX00079', 'NIX00080', 'NIX00081']
}, {
    'Attenuator': ['NIXD0050', 'NIXD0051', 'NIXD0088', 'NIXD0089']
}, {
    'Attenuator': ['NIXD0068', 'NIXD0069', 'NIX00076', 'NIX00094', 'NIXD0075']
}, {
    'HV and depolarization': ['NIXD0023', 'NIXD0090', 'NIXD0091']
}, {
    'HV and depolarization': ['NIXD0023', 'NIXD0066', 'NIXD0067', 'NIXD0074']
}, {
    'High voltage': ['NIXD0052', 'NIXD0053']
}, {
    'Detector temperature sensors':
    ['NIXD0054', 'NIXD0055', 'NIXD0056', 'NIXD0057']
}, {
    'Detector': [
        'NIXD0082',
        'NIXD0083',
        'NIXD0084',
        'NIXD0085',
        'NIXD0070',
        'NIXD0058',
    ]
}, {
    'Trigger rates': ['NIXD0023', 'NIX00072', 'NIX00073']
}]

#delta_time between packets

figsize = (10, 7)


def get_delta_time(time_array):
    last_ts = time_array[0]
    delta_ts = []
    for t in time_array[1:]:
        delta_ts.append(t - last_ts)
        last_ts = t
    return delta_ts


class Plugin(object):
    def __init__(self, packets=[]):
        self.packets = packets

    def run(self, pdf):
        print('Number of packets : {}'.format(len(self.packets)))

        param_values = sdt.Packet.merge(self.packets, SPIDs)
        #with PdfPages(filename) as pdf:
        fig = None
        fig = plt.figure(figsize=figsize)
        plt.axis('off')
        try:
            title = ' Housekeeping data SPID(s): {} ,start: {} stop: {}'.format(
                str(SPIDs), param_values['UTC'][0], param_values['UTC'][-1])
        except KeyError:
            title = ' Housekeeping data SPID(s): {} ,start: {} stop: {}'.format(
                str(SPIDs), param_values['unix_time'][0],
                param_values['unix_time'][-1])
        plt.text(0.5, 0.5, title, ha='center', va='center')

        pdf.savefig()
        plt.close()
        fig.clf()

        fig = plt.figure(figsize=figsize)
        title = "Timestamps"
        fig.suptitle(title, fontsize=14, fontweight='bold')
        ax = fig.add_subplot(2, 1, 1)
        ax.plot(param_values['unix_time'])
        ax.set_xlabel('Packet #')
        ax.set_ylabel('Time (s)')
        ax.set_title('Packet time')

        ax = fig.add_subplot(2, 1, 2)
        ax.plot(get_delta_time(param_values['unix_time']))
        ax.set_xlabel('Packet #')
        ax.set_ylabel('Delta T(s)')
        ax.set_title('Packet timestamp difference')
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        pdf.savefig()
        plt.close()
        fig.clf()
        for group in groups:
            fig = plt.figure(figsize=figsize)
            title = list(group.keys())[0]
            group_parameters = list(group.values())[0]
            num_plots = len(group_parameters)
            nrows = 0
            ncols = 1
            if num_plots <= 3:
                nrows = num_plots
            else:
                nrows = int(math.sqrt(num_plots))
            ncols = math.ceil(num_plots / nrows)
            fig.suptitle(title, fontsize=14, fontweight='bold')
            for ifig, pname in enumerate(group_parameters):
                ax = fig.add_subplot(nrows, ncols, ifig + 1)
                text_cal = STIX_IDB.get_textual_mapping(pname)
                value = param_values[pname]
                if len(param_values['unix_time']) == len(value):
                    ax.step(param_values['unix_time'], value, where='mid')
                else:
                    ax.plot(value)
                ax.tick_params(axis='x', which='major', pad=10)
                fig.tight_layout()
                ax.set_xlabel('Time (s)')
                desc = stix_desc.get_parameter_desc(pname)
                ax.set_title('{} ({})'.format(desc, pname))
                if text_cal:
                    if len(text_cal[0]) < 8:
                        ax.set_yticks(text_cal[0])
                        ax.set_yticklabels(text_cal[1])

                else:
                    ax.set_ylabel('value')
            fig.tight_layout()
            fig.subplots_adjust(top=0.9)
            pdf.savefig()
            plt.close()
            fig.clf()
