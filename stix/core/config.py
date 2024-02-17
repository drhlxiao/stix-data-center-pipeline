#!/usr/bin/python3
"""
    stix parser configuration manager
    created on Oct. 20, 2020
"""
import os
import json
import numpy as np
from dateutil import parser as dtparser
from stix.core import logger
logger = logger.get_logger()

ASW_VERSION = 186
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
                "aspect_l2_path_pattern":"/data/pub099/fits/L2/*/*/*/*/*-aux-*.fits",
                "flare_lc_snapshot_path": "/data/flare_lc",
                "calibration_report_path": "/data/calibration/",
                "level1_products_path": "/data/level1/",
                "level2_products_path": "/data/level2/",
                "ngnix_cache": "/data/nginx/cache/*",
                "goes_lc_path": "/data/goes/",
                "flare_pipeline_path":"/data/flare_pipeline"
                },
            "asw_version": ASW_VERSION
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
                "validityPeriod": ["2020-12-28T00:00:00", "2021-12-10T14:00:00"]
                },
            "183":{
                "filename": "/data/pub/data/idb/idb_v2.26.36.sqlite",
                "version": "2.26.36",
                "validityPeriod": ["2021-12-10T14:00:00", "2024-02-05:00:00"]
                },
            "186":{
                "filename": "/data/pub/data/idb/idb_v2.26.37.sqlite",
                "version": "2.26.37",
                "validityPeriod": ["2024-02-05T00:00:00", "2030-02-07:00:00"]
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
instrument_config={
        'default_ebins':[0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70,
                      76, 84, 100, 120, 150, np.inf],
        'ebins_history': [
        {
            'time_range': [1690761600, 1699140003],  # 2023-07-31T00:00:00 to 2023-11-05
            'energy_bins': [0.0, 4.0, 5.5, 6.3, 6.56, 6.8, 7.3, 8.0, 9.0,
                            10.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0,
                            28.0, 32.0, 36.0, 40.0, 45.0, 50.0, 56.0, 63.0, 70.0, 76.0, 84.0, 100.0, 120.0, 150.0, np.inf]
        },
        {
            'time_range': [1675036800, 1675353600],  # 2023-01-30 to 2023-02-02T16:00
            'energy_bins': [0, 4, 4.45, 4.95, 5.45, 5.95, 6.45, 6.95, 7.35, 7.75, 8.25, 8.75, 9.25, 10,
                            10.5, 11, 12, 13, 15, 18, 21, 25, 28, 32, 36, 43, 50, 59, 70, 84, 110, 150, np.inf]
        }
    ],
        'ql_trig_default_scale_factor':30,
        'sci_trig_default_scale_factor':30,

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
