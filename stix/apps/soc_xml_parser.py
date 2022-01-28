#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import (absolute_import, unicode_literals)
import argparse
import pprint
import binascii
from io import BytesIO
import xmltodict
from core import idb
from core import global
from core import writer_sqlite as stw
from core import logger
from core import telemetry_parser as tm_parser

LOGGER = logger.LOGGER


def parse_esa_xml_file(in_filename, out_filename=None, selected_spid=0):
    packets = []
    with open(in_filename) as fd:
        doc = xmltodict.parse(fd.read())
        for e in doc['ns2:ResponsePart']['Response']['PktRawResponse'][
                'PktRawResponseElement']:
            packet = {'id': e['@packetID'], 'raw': e['Packet']}
            packets.append(packet)
        num_packets = 0
        num_fix_packets = 0
        num_variable_packets = 0
        num_bytes_read = 0
        st_writer = stw.writer(out_filename)
        st_writer.register_run(in_filename)
        total_packets = 0
        for packet in packets:
            data_hex = packet['raw']
            data_binary = binascii.unhexlify(data_hex)
            in_file = BytesIO(data_binary[76:])
            status, header, parameters, param_type, param_desc, num_bytes_read = tm_parser.parse_one_packet(
                in_file, LOGGER)
            total_packets += 1
            LOGGER.pprint(header, parameters)
            if param_type == 1:
                num_fix_packets += 1
            elif param_type == 2:
                num_variable_packets += 1
            if header and parameters:
                st_writer.write_parameters(parameters)

        LOGGER.info('{} packets found in the file: {}'.format(
            total_packets, in_filename))
        LOGGER.info('{} ({} fixed and {} variable) packets processed.'.format(num_packets,\
                num_fix_packets,num_variable_packets))
        LOGGER.info('Writing parameters to file {} ...'.format(out_filename))
        st_writer.done()
        LOGGER.info('Done.')


def main():
    in_filename = 'test/stix.xml'
    out_filename = 'out.db'
    sel_spid = 0
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--in", required=True, help="input file")
    ap.add_argument("-o", "--out", required=False, help="output file")
    ap.add_argument(
        "-s",
        "--sel",
        required=False,
        help="only select packets of the given SPID")

    args = vars(ap.parse_args())
    if args['out'] is not None:
        out_filename = args['out']
    if args['sel'] is not None:
        sel_spid = int(args['sel'])
    in_filename = args['in']
    LOGGER.info('Input file', in_filename)
    LOGGER.info('Output file', out_filename)
    parse_esa_xml_file(in_filename, out_filename, sel_spid)


if __name__ == '__main__':
    main()
