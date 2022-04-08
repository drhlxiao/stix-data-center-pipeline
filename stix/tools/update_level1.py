from astropy.io import fits
from stix.spice import time_utils as sdt
from stix.spice import spice_manager as spm
from astropy.time import Time
import glob
import json
def update(fname):
   with open(fname, "r+") as jsonFile:
       data = json.load(jsonFile)
       data["status"] = "OK"
       jsonFile.seek(0)  # rewind
       json.dump(data, jsonFile)
       jsonFile.truncate()
for fname in glob.glob('/data/level1/*.json'):
    print(fname)
    update(fname)
