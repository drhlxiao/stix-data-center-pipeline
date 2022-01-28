#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @title        : xml_to_bin.py
# @description  : read packets from an ESA xml file and write them to STIX raw data format
# @date         : March. 28, 2019
from __future__ import (absolute_import, unicode_literals)
import sys
import binascii
import xmltodict


def convert_xml_to_bin(in_filename, out_filename):
    with open(in_filename) as fin, \
            open(out_filename,'wb') as fout :
        doc = xmltodict.parse(fin.read())
        num_packet = 0
        for packet in doc['ns2:ResponsePart']['Response']['PktRawResponse'][
                'PktRawResponseElement']:
            data_hex = packet['Packet']
            data_binary = binascii.unhexlify(data_hex)
            fout.write(data_binary[76:])
            num_packet += 1
        print(
            ' %d packets have been written to %s' % (num_packet, out_filename))


def main():
    if len(sys.argv) != 3:
        print('xml_to_bin  <INPUT> <OUTPUT>')
    else:
        convert_xml_to_bin(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()
