import sys
sys.path.append('.')
import requests
from datetime import datetime
from stix.core import mongo_db as db
from stix.spice import stix_datetime
from stix.core import config
from stix.wiki import plot_orbit,aia, stix_lightcurves, goes
mdb = db.MongoDB()

def process_flares_for_file(file_id):
    fdb=mdb.get_collection('flares_tbc')
    flares=fdb.find({'run_id':file_id, 'hidden':False})
    folder=config.get_config('joint_obs')
    if not flares:
        print(f'Flare {file_id} not file in db!')
        return
    for doc in flares:
        start_unix=doc['start_unix']
        end_unix=doc['end_unix']
        peak_utc=doc['peak_utc']
        start_utc=stix_datetime.unix2utc(start_unix)
        end_utc=stix_datetime.unix2utc(end_unix)
        flare_id=doc['flare_id']
        _id=doc['flare_id']
        print('generate location')
        emph=plot_orbit.plot_solo_location(folder,_id, flare_id,peak_utc,  overwrite=False)
        print('processing aia')
        aia.plot_aia(folder,_id, flare_id, peak_utc,  wavelen=131, overwrite=False)
        print('processing goes')
        aia.plot_aia(folder,_id, flare_id, peak_utc,  wavelen=131, overwrite=False)
        goes.plot_goes(folder,_id, flare_id, start_utc, end_utc, overwrite=False)
        try:
            light_time=emph['light_time_diff'][0]
        except:
            light_time=0
        stix_lightcurves.plot_stix_lc(folder, _id, flare_id, start_utc, end_utc,  overwrite=False, T0_utc=peak_utc, light_time=light_time)

if __name__=='__main__':
    if len(sys.argv)==2:
        file_id=int(sys.argv[1])
        process_flares_for_file(file_id)
    else:
        print('file id not provided')
