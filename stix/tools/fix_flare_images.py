
from stix.core import mongo_db as db
import astropy.units as u
from datetime import datetime
from stix.spice import solo
from stix.spice import time_utils as ut
mdb = db.MongoDB()
flare_images_db= mdb.get_collection('flare_images')


for doc in flare_images_db.find():
    obs_utc=doc['utc_range'][0]

    print(doc['_id'])
    try:
        stix_aux=solo.SoloEphemeris.get_solar_limb_stix_fov(obs_utc)
        sun_center=stix_aux['sun_center']
        pointing_valid=True
    except :
        sun_center=[0,0]
        pointing_valid=False
    start=ut.utc2filename(obs_utc)
    rsun=np.degrees(np.arctan(stix_aux['rsun']/stix_aux['solo_sun_r']))

    flare_images_db.update_one({'_id':doc['_id']},
            {'$set':
    {'idl_status':'',
        'num_idl_calls':0,
        'fits':{},
        'figs':{},
        'idl_config':{
            'folder':f'/data/quicklook/flare_images/{doc["unique_id"]}',
            'prefix':f'stix_image_sci_{doc["bsd_id"]}_uid_{doc["unique_id"]}_{doc["energy_range"][0]}_{doc["energy_range"][1]}_{start}_{doc["_id"]}',
            'fwdfit_shape': 'ellipse' if doc['energy_range'][1]<15 else 'multi'
            },
        "creation_time":datetime.now(),
        "aux.dsun":stix_aux['solo_sun_r'],
        "aux.rsun":rsun*3600,
        "aux.sun_center": sun_center,
        "aux.has_pointing": pointing_valid
        }})

