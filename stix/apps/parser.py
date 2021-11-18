#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : parser.py
# @description  : STIX TCTM packet parser
# @author       : Hualin Xiao
# @date         : Feb. 11, 2019
#

import argparse
import os
import sys
from stix.core import config
from stix.core import stix_logger, stix_idb, stix_parser
logger = stix_logger.get_logger()
MONGODB_CONFIG=config.get_config('pipeline.mongodb')


def main():

    ap = argparse.ArgumentParser()
    required = ap.add_argument_group('Required arguments')
    optional = ap.add_argument_group('Optional arguments')


    required.add_argument(
        "-i",
        dest='input',
        required=True,
        nargs='?',
        help="Input raw data filename.")
    optional.add_argument(
        "-o",
        dest='output',
        default=None,
        required=False,
        help="Output filename. ")

    optional.add_argument(
        "--idb",
        dest='idb',
        default=None,
        required=False,
        help="IDB sqlite3 filename. ")


    optional.add_argument(
        "-t",
        dest='input_type',
        default=None,
        choices=('bin', 'ascii', 'xml', 'hex'),
        help=
        "Input file type. Four types (bin, hex, ascii or xml) are supported.")

    optional.add_argument(
        "--wdb",
        dest='wdb',
        default=False,
        action='store_true',
        help='Write decoded packets to local MongoDB.')

    optional.add_argument(
        "-s",
        dest='source',
        default='',
        required=False,
        help='data source: GU or PFM')


    optional.add_argument(
        "--db-host",
        dest='db_host',
        default='localhost',
        help='MongoDB host IP.')
    optional.add_argument(
        "--db-port",
        dest='db_port',
        default=MONGODB_CONFIG['port'],
        type=str,
        help='MongoDB host port.')

    optional.add_argument(
        "--db-user", dest='db_user', 
        default=MONGODB_CONFIG['user'], help='MongoDB username.')
    optional.add_argument(
        "--db-pwd", dest='db_pwd', default=MONGODB_CONFIG['password'], help='MongoDB password.')
    optional.add_argument(
        "-m", default='', dest='comment', required=False, help="comment")

    optional.add_argument(
        '--SPID',
        nargs='*',
        dest='SPID',
        action="store",
        default=[],
        type=int,
        help='Only to parse the packets of the given SPIDs.')

    optional.add_argument(
        '--services',
        nargs='*',
        dest='services',
        action="store",
        default=[],
        type=int,
        help='Only to parse the packets of the given service types.')
    optional.add_argument(
        "-v",
        dest="verbose",
        default=5,
        required=False,
        help="Logger verbose level",
        type=int)
    optional.add_argument(
        "-l",
        "--log",
        dest='logfile',
        default=None,
        required=False,
        help="Log filename")

    optional.add_argument(
        "--no-S20",
        dest='S20_excluded',
        action='store_true',
        default=True,
        help="to exclude S20 packets")

    args = vars(ap.parse_args())


    logger.set_logger(args['logfile'], args['verbose'])

    if not os.path.exists(args['input']):
        logger.error("File {} doesn't exist.".format(args['input']))
        return

    if args['idb']:
        idb_instance = stix_idb.stix_idb(args['idb'])
    selected_spids = args['SPID']
    selected_services = args['services']
    comment=args['comment']
    process_single_file(args['input'], args['input_type'],args['output'], args['SPID'],
            args['services'],args['comment'], args['S20_excluded'],args['wdb'], args['db_host'],
            args['db_port'], args['db_user'], args['db_pwd'], args['source'])

def process_single_file(filename, filetype, output_filename, selected_spids,selected_services ,comment, S20_excluded,
        wdb, db_host,db_port, db_user, db_pwd, source):
    parser = stix_parser.StixTCTMParser()
    parser.set_packet_filter(selected_services, selected_spids)
    if S20_excluded:
        parser.exclude_S20()
    if output_filename:
        parser.set_packet_buffer_enabled(False)
        parser.set_store_binary_enabled(False)
        #parser.set_pickle_writer(output_filename, comment)
    if wdb:
        parser.set_store_binary_enabled(False)
        parser.set_packet_buffer_enabled(False)
        parser.config_mongodb(db_host,db_port,db_user, db_pwd,comment, filename)
    parser.parse_file(filename, filetype) 
    parser.done()


if __name__ == '__main__':
    main()
