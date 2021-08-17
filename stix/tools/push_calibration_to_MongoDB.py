import os
import sys
import argparse
from matplotlib.backends.backend_pdf import PdfPages

from core import stix_logger
from core import stix_parser
from core import stix_idb
from core import stix_calibration

from utils import mongo_db
STIX_LOGGER = stix_logger.stix_logger()
STIX_MDB= mongo_db.MongoDB()


def scan_all():
    runs=STIX_MDB.select_all_runs(1)
    collection=STIX_MDB.get_collection_calibration()
    calibration=stix_calibration.StixCalibration(collection)
    for run in runs:
        packets=STIX_MDB.select_packets_by_run(run['_id'])
        print(run['_id'],len(packets))
        for packet in packets:
            calibration.capture(run['_id'],packet['_id'],packet)

scan_all()





