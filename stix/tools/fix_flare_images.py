
from stix.core import mongo_db as db
import astropy.units as u
from datetime import datetime
from stix.spice import solo
mdb = db.MongoDB()
flare_images_db= mdb.get_collection('flare_images')


for doc in flare_images_db.find():
    obs_utc=doc['utc_range'][0]
    try:
        stix_aux=solo.SoloEphemeris.get_solar_limb_stix_fov(obs_utc)
        sun_center=stix_aux['sun_center']
        pointing_valid=True
    except :
        sun_center=[0,0]
        pointing_valid=False
    flare_images_db.update_one({'_id':doc['_id']},
            {'$set':
    {'idl_status':'',
        'num_idl_calls':0,
        'fits':{},
        'figs':{},
        'idl_config':{
            'folder':f'/data/quicklook/flare_images/{doc["unique_id"]}',
            'prefix':f'stix_image_sci_{doc["bsd_id"]}_uid_{doc["unique_id"]}_{doc["energy_range"][0]}_{doc["energy_range"][1]}_{doc["utc_range"][0]}_{doc["_id"]}',
            'fwfit_shape': 'ellipse' if doc['energy_range'][1]<15 else 'multi'
            },
        "creation_time":datetime.now(),
        "aux.dsun":(doc["aux"]['solo_sun_r']*u.au).to(u.m).value,
        "aux.sun_center": sun_center,
        "aux.has_pointing": pointing_valid
        }})

