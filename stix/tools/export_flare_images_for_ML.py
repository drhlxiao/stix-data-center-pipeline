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
def get_data(filename):
    m=sunpy.map.Map(filename)
    xyc=ndimage.center_of_mass(m.data)
    sc=m.pixel_to_world(xyc[0]*u.pix,xyc[1]*u.pix)
    center=[sc.Tx.to(u.arcsec).value, sc.Ty.to(u.arcsec).value]

    return m.data, center


for doc in flare_images_db.find({'total_counts':{'$gt':5000}, 'duration':{'$lt':900}, 
    'run_type':'auto'}):
    pixel_counts=np.array(doc['pixel_counts'])

    try:
        fn=doc['fits']['image_em']
        em, em_flare_center=get_data(fn)
        if np.sum(em)>0:
            em_fname=f'/home/xiaohl/flare_images/em_{doc["_id"]}.npz'
            np.savez(em_fname,  energy_range=np.array(doc['energy_range']), image=em.astype(np.float32),
                    counts=pixel_counts,  center= np.array(em_flare_center, dtype=np.float32))
    except:
        pass


    try:
        fn=doc['fits']['image_clean']
        clean, clean_flare_center=get_data(fn)
        if np.sum(clean)>0:
            clean_fname=f'/home/xiaohl/flare_images/clean_{doc["_id"]}.npz'
            np.savez(clean_fname,  energy_range=np.array(doc['energy_range']), image=clean.astype(np.float32),
                    counts=pixel_counts,  center= np.array(clean_flare_center, dtype=np.float32))
            print(clean_fname)

    except:
        pass


