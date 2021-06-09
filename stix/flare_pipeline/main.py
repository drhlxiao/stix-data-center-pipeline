import sys
sys.path.append('.')
import requests
from datetime import datetime
from stix.core import mongo_db as db
from stix.spice import stix_datetime
from stix.core import config
from stix.flare_pipeline import plot_orbit,aia,solo_eui, stix_lightcurves, goes, wiki_creator, aia_solo_view
mdb = db.MongoDB()
margin=1200
threshold=600
OVERWRITE=False

def to_earth_utc(utc, tdiff):
    return stix_datetime.unix2utc(stix_datetime.utc2unix(utc)+tdiff)

def process_flares_for_file(file_id, overwrite=OVERWRITE):
    print(f'processing for file {file_id}')
    fdb=mdb.get_collection('flares')
    flares=fdb.find({'run_id':file_id, 'hidden':False})
    folder=config.get_config('joint_obs')
    if not flares:
        print(f'Flare not found in  {file_id} !')
        return
    for doc in flares:
        start_unix=doc['start_unix']-margin
        end_unix=doc['end_unix']+margin
        peak_utc=doc['peak_utc']
        peak_counts=doc['peak_counts']
        if peak_counts<=threshold:
            print(f"Ignored flares {doc['_id']}, peak counts < {threshold}")
            continue
        start_utc=stix_datetime.unix2utc(start_unix)
        end_utc=stix_datetime.unix2utc(end_unix)
        flare_id=doc['flare_id']
        _id=doc['_id']
        print(_id)
        print('generate location')
        eph=plot_orbit.plot_solo_location(folder,_id, flare_id,peak_utc,  overwrite=overwrite)
        if eph:
            mdb.update_flare_field(_id, 'ephemeris',
                    {
                        'x':eph['x'][0],
                        'y':eph['y'][0],
                        'z':eph['z'][0],
                        'lt_diff':eph['light_time_diff'][0],
                        'earth_sun_solo_angle':eph['earth_sun_solo_angle'][0],
                        'sun_solo_r':eph['sun_solo_r'][0]
                })
        try:
            light_time=eph['light_time_diff'][0]
        except Exception as e:
            print(e)
            light_time=0
        stix_lightcurves.plot_stix_lc(folder, _id, flare_id, start_utc, end_utc,  overwrite=overwrite, peak_utc=peak_utc, light_time=light_time)

        earth_start_utc=to_earth_utc(start_utc,light_time) 
        earth_end_utc=to_earth_utc(end_utc,light_time) 
        earth_peak_utc=to_earth_utc(peak_utc,light_time) 
        try:
            goes.plot_goes(folder,_id, flare_id,earth_start_utc,earth_end_utc, earth_peak_utc, overwrite=overwrite)
        except Exception as e:
            raise
            print("error when creating goes LC")
            print(e)

        try:
            aia.plot_aia(folder,_id, flare_id, earth_peak_utc,  131, overwrite)
        except Exception as e:
            print("error when creating AIA 131 image")
            print(e)
        try:
            aia_solo_view.plot(folder,_id, flare_id, earth_peak_utc,  1600, overwrite)
        except Exception as e:
            print("error when creating AIA 1600 image")
            print(e)
        try:
            sun_angular_diameter_arcmin=eph.get('sun_angular_diameter_arcmin',60)
            solo_eui.plot_eui(folder,_id, flare_id, start_utc, end_utc, sun_angular_diameter_arcmin, overwrite)
        except Exception as e:
            print("error when creating EUI image")
            print(e)
        try:
            wiki_creator.wiki_bot.touch_wiki_for_flare(_id)
        except Exception as e:
            print("error when creating wiki")
            print(e)
            

if __name__=='__main__':
    if len(sys.argv)==2:
        file_id=int(sys.argv[1])
        process_flares_for_file(file_id)
    elif len(sys.argv)==3:
        file_id=int(sys.argv[1])
        file_id_end=int(sys.argv[2])
        for i in range(file_id, file_id_end+1):
            process_flares_for_file(i)
    else:
        print('file id not provided')
            

