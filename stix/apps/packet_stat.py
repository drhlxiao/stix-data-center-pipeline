#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : packet_stat.py
# @description  : packet statistics
# @author       : Hualin Xiao
# @date         : Feb. 11, 2019
#

from __future__ import (absolute_import, unicode_literals)
import argparse
import pprint
from core import idb
from core import global
from core import writer
from core import logger
from core import parser
from core import plotter

LOGGER = logger.LOGGER
_idb = idb._idb


def get_packet_type_stat_text(spid_list):
    counter = [{
        'spid': x,
        'counts': spid_list.count(x)
    } for x in set(spid_list)]
    text = '#    * SPID  * Service Type * Service SubType *Number of packets*  Description\n'

    sorted_counter = sorted(counter, key=lambda k: k['spid'])
    descr = dict()

    for i, item in enumerate(sorted_counter):
        spid = item['spid']
        counts = item['counts']
        desc = ''
        pid_type = 0
        pid_stype = 0
        try:
            row = _idb.get_spid_info(spid)
            pid_type = row[0][1]
            pid_stype = row[0][2]
            desc = str(row[0][0])
            descr[spid] = desc
        except:
            pass
        text += (
            '{:<3}  * {:<6} *  {:<10} * {:<14}  *  {:<16}  *  {:<32}\n'.format(
                i, spid, pid_type, pid_stype, counts, str(desc)))
    text += '-' * 70
    return text


def analyze_raw_file(in_filename, noplot=False):
    """
    Parse stix raw packets 
    Args:
     in_filename: input filename
     out_filename: output filename
     selected_spid: filter data packets by  SPID. 0  means to select all packets
    Returns:
    """
    spid_list = []
    pid = 0
    with open(in_filename, 'rb') as in_file:
        num_packets = 0
        num_bytes_read = 0
        timestamps = []
        while True:
            status, header, header_raw, application_data_raw, bytes_read = parser.read_one_packet(
                in_file, LOGGER)
            num_bytes_read += bytes_read
            if status == global.NEXT_PACKET:
                continue
            if status == global.EOF:
                break

            spid = header['SPID']
            spid_list.append(spid)
            timestamps.append(header['time'])

        print('-' * 70)
        print('Data file       : {}'.format(in_filename))
        print('Total bytes read: {} kB'.format(num_bytes_read / 1024))
        print('Total packets   : {}'.format(len(spid_list)))
        text = get_packet_type_stat_text(spid_list)
        print(text)

        if not noplot:
            print('preparing timeline plot ...')
            plotter.plot_packet_header_timeline(timestamps, spid_list)


def main():
    in_filename = 'test/stix.dat'
    noplot = False
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--in", required=True, help="input file")
    ap.add_argument("-o", "--out", required=False, help="output file")
    ap.add_argument(
        '-np', action='store_true', help='Don\'t create a timeline plot')

    args = vars(ap.parse_args())
    if args['np'] is not None:
        noplot = args['np']
    in_filename = args['in']
    analyze_raw_file(in_filename, noplot)


if __name__ == '__main__':
    main()
