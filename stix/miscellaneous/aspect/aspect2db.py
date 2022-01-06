# read aspect fits file and push data to files
import os
import numpy
from astropy.io import fits
from stix.spice import stix_datetime 
#from stix.core import config
#level2_products_path = config.get_config(
#    'pipeline.daemon.level2_products_path')

mdb = db.MongoDB()
def create_fits_doc(filename, path, start_unix, end_unix, creation_utc)
    fits_id=mdb.get_next_fits_id()
    filename=os.path.basename(filename)
    doc={"_id" : fits_id,
    "packet_id_start" : -1,
    "packet_id_end" : -1,
    "packet_spid" : -1,
    "num_packets" : -1 ,
    "file_id" : -1,
    "product_type" : "science",
    "product_group" : "aspect",
    "complete" : true,
    "version" : 1,
    "level" : "L2",
    "creation_time" : stix_datetime.utc2datetime(creation_utc), 
    "path" : path,
    "data_start_unix" : start_unix,
    "data_end_unix" : end_unix,
    "filename" : filename,
    "file_size" : 28800
    }
    mdb.write_fits_index_info(doc)

def monitor_directory(folder='.'):
    pass


def process(folder, filename):
    hdul = fits.open(fname):
    asp_fits_db=mdb.get_collection('aspect_fits')
    asp_solution_db=mdb.get_collection('aspect_solutions')


    header=hdul['Primary'].header
    metadata={
            'creation_date':header['DATE'],
            'data_begin': header['DATE_BEG'],
            'data_end': header['DATE_END'],
            'start_unix': stix_datetime.utc2unix(header['DATE_BEG']),
            'end_unix': stix_datetime.utc2unix(header['DATE_END']),
            'filename':os.path.filename
            }
    ['DATE','DATE_BEG','DATE_END']
    for row in hdul['DATA'].data.tolist():
