#export flare images for machine learning
import numpy as np
from stix.core import mongo_db as db
import astropy.units as u
from astropy.io import fits 
from datetime import datetime
from stix.spice import solo
from stix.spice import time_utils as ut
from scipy import ndimage
import pickle
import h5py
import sunpy
import sunpy.map
mdb = db.MongoDB()
flare_images_db= mdb.get_collection('flare_images')
docs=list(flare_images_db.find({'ospex_meta':True, 'energy_range.0':4, 'energy_range.1':10, 'flare_aux':{'$exists':True}}))
with open('flare_image_db.pkl') as f:
    pickle.dump(docs, f)

