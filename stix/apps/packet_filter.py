#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : packet_filter.py
# @description  : a script to filter stix TM packets by SPID
# @author       : Hualin Xiao
# @date         : March 12, 2019
#

from __future__ import (absolute_import, unicode_literals)
import argparse
import pprint
from core import idb
from core import global
from core import parser

_telmetry_parser = parser.StixTCTMParser()


def packet_filter(in_filename, out_filename, selected_spid):
    """
    Parse stix raw packets 
    Args:
     in_filename: input filename
     out_filename: output filename
     selected_spid: filter data packets by  SPID. 0  means to select all packets
    Returns:
    """
    with open(in_filename, 'rb') as fin, \
            open(out_filename, 'wb') as fout:
        num_selected_packets = 0
        num_total_packets = 0
        while True:
            result = _telmetry_parser.read_packet(fin)
            status = result['status']
            header = result['header']
            header_raw = result['header_raw']
            app_raw = result['app_raw']
            num_bytes_read = result['num_read']
            if status == global._next_packet:
                continue
            if status == global._eof:
                break
            spid = header['SPID']
            num_total_packets += 1
            if spid == selected_spid:
                num_selected_packets += 1
                fout.write(header_raw)
                fout.write(app_raw)


def main():
    in_filename = 'test/stix.dat'
    out_filename = 'export.dat'
    selected_spid = 0
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--in", required=True, help="input file")
    ap.add_argument("-o", "--out", required=True, help="output file")
    ap.add_argument(
        "-s", "--sel", required=True, help="select packets by SPID")

    args = vars(ap.parse_args())
    if args['out'] is not None:
        out_filename = args['out']
    if args['sel'] is not None:
        selected_spid = int(args['sel'])
    in_filename = args['in']
    packet_filter(in_filename, out_filename, selected_spid)


if __name__ == '__main__':
    main()
