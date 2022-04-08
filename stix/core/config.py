#!/usr/bin/python3
"""
    stix parser configuration manager
    created on Oct. 20, 2020
"""
import os
import json
from dateutil import parser as dtparser
from stix.core import logger
logger = logger.get_logger()

ASW_VERSION = 183
HTTP_PREFIX = 'https://pub023.cs.technik.fhnw.ch'
parser_config = {
        "pipeline": {
            "mongodb": {
                "host": "localhost",
                "user": "",
                "password": "",
                "port": 27017
                },
            "daemon": {
                "data_source": {
                    "GU": ["/home/xiaohl/data/*.ascii"],
                    "PFM": [
                        "/data/gfts/solmoc/from_moc/*.xml",
                        "/data/gfts/solmoc/from_edds/tm/*.xml",
                        "/data/gfts/solmoc/from_edds/tc/*.xml",
                        "/data/gfts/solmoc/from_moc/*ascii"
                        ]
                    },
                "log_path": "/data/log/",
                "flare_images": "/data/quicklook/flare_images",
                "notification": {
                    "file": "/data/log/message.log"
                    },
                "fits_path": "/data/fits",
                "flare_lc_snapshot_path": "/data/flare_lc",
                "calibration_report_path": "/data/calibration/",
                "level1_products_path": "/data/level1/",
                "level2_products_path": "/data/level2/",
                "ngnix_cache": "/data/nginx/cache/*",
                "goes_lc_path": "/data/goes/",
                "flare_pipeline_path":"/data/flare_pipeline"
                },
            "asw_version": 183
            },
        "ASW": {
            "179": {
                "filename": "/data/pub/data/idb/idb.sqlite",
                "version": "2.26.33",
                "aswVersion": 179,
                },
            "181": {
                "filename": "/data/pub/data/idb/idb_v2.26.35.sqlite",
                "version": "2.26.34",
                "validityPeriod": ["2020-12-28T00:00:00", "2020-12-28T00:00:00"]
                },
            "183":{
                "filename": "/data/pub/data/idb/idb_v2.26.36.sqlite",
                "version": "2.26.36",
                "validityPeriod": ["2020-12-28T00:00:00", "2020-12-28T00:00:00"]
                }
            },
        'joint_obs':     '/data/flares/',
"spice": "/data/pub/data/spice/latest/kernels",
        'mailer':    {
                'sender':  'noreply@fhnw.ch',
                'user': '',
                'pwd': '',
                'server': 'lmailer.fhnw.ch',
                'port': 465
                }



}


def config(key):
    return parser_config.get(key, '')


def get_config(key=None):
    # get configuration value
    # For example:  get_config(pipeline.daemon.fits_path)
    if not key:
        return parser_config

    if '.' in key:
        result = parser_config
        try:
            for item in key.split('.'):
                result = result[item]
            return result
        except IndexError or ValueError:
            logger.error(f'Can not find  {key} in config')
            return None
    return parser_config.get(key, '')


def get_idb(asw_version=None):
    if not asw_version:
        asw_version = ASW_VERSION

    fname=parser_config["ASW"][str(asw_version)]["filename"]
    #print(fname)
    return fname


def get_spice_folder():
    return parser_config['spice']


#print(config)
#print(get_config('pipeline.mongodb.host'))
#print(get_idb('2020-10-01T00:00:00'))
#print(get_spice('2020-10-01T00:00:00'))
#pprint(config)
