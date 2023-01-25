import sys
import os
import argparse
from collections import defaultdict
from datetime import datetime, timedelta
from itertools import chain
from pathlib import Path

from stix.core import datatypes as sdt
from stix.fits.io.processors import FitsL1Processor
from stix.fits.io import hk_fits_writer as hkw
from stix.spice import time_utils as st
from stix.fits.products.housekeeping import MiniReport, MaxiReport
from stix.fits.products.quicklook import LightCurve, Background, Spectra, Variance, \
    FlareFlagAndLocation, CalibrationSpectra, TMManagementAndFlareList

from stix.fits import fits_creator
from stix.core import mongo_db, logger

logger = logger.get_logger()

db = mongo_db.MongoDB()

bsd_db = db.get_collection('bsd')
fits_db = db.get_collection('fits')
packets_db = db.get_collection('packets')

docs = bsd_db.find({'SPID':54143}).sort('_id',1)


for doc in docs:
    #unique_id=doc['unique_id']
    fits_creator.purge_fits_for_science_request(unique_id)




