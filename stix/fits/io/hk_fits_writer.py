import os
from datetime import datetime
from stix.spice.datetime import datetime_to_scet
def generate_primary_header(filename, scet_coarse, scet_fine, obs_beg, obs_avg, obs_end, version):
    """
    Generate primary header cards.

    Parameters
    ----------
    filename : str
        Filename
    scet_coarse : int
        SCET coarse time
    scet_fine : int
        SCET fine time
    obs_beg : datetime.datetime
        Begging of observation
    obs_avg : datetime.datetime
       Averagea of observation
    obs_end : datetime.datetime
        End of observation

    Returns
    -------
    tuple
        List of header cards as tuples (name, value, comment)
    """
    headers = (
        # Name, Value, Comment
        ('TELESCOP', 'SOLO/STIX', 'Telescope/Sensor name'),
        ('INSTRUME', 'STIX', 'Instrument name'),
        ('OBSRVTRY', 'Solar Orbiter', 'Satellite name'),
        ('FILENAME', filename, 'FITS filename'),
        ('DATE', datetime.now().isoformat(timespec='milliseconds'),
         'FITS file creation date in UTC'),
        ('OBT_BEG', f'{scet_coarse[0]}:{scet_fine[0]}'),
        ('OBT_END', datetime_to_scet(obs_end)),
        ('TIMESYS', 'UTC', 'System used for time keywords'),
        ('LEVEL', 'L1A', 'Processing level of the data'),
        ('ORIGIN', 'STIX Team, FHNW', 'Location where file has been generated'),
        ('CREATOR', 'STIX-SWF', 'FITS creation software'),
        ('VERSION', version, 'Version of data product'),
        ('OBS_MODE', 'Nominal '),
        ('VERS_SW', 1, 'Software version'),
        ('DATE_OBS', obs_beg.isoformat(timespec='milliseconds'),
         'Start of acquisition time in UT'),
        ('DATE_BEG', obs_beg.isoformat(timespec='milliseconds')),
        ('DATE_AVG', obs_avg.isoformat(timespec='milliseconds')),
        ('DATE_END', obs_end.isoformat(timespec='milliseconds')),
        #('OBS_TYPE', 'LC'),
        # TODO figure out where this info will come from
        ('SOOP_TYP', 'SOOP'),
        ('OBS_ID', 'obs_id'),
        ('TARGET', 'Sun')
    )
    return headers

def generate_filename(product,  unique_id, product_type,  version, run_type):
    """
    Generate fits file name with SOLO conventions.

    Parameters
    ----------
    level : str
        Data level e.g L0, L1, L2
    product_name : str
        Name of the product eg. lightcruve spectra
    observation_date : datetime.datetime
        Date of the observation
    version : int
        Version of this product

    Returns
    -------
    str
        The filename
    """
    #dateobs = observation_date.strftime("%Y%m%dT%H%M%S")
    #return f'solo_{level}_stix-{product_name.replace("_", "-")}_{dateobs}_{unique_id:05d}.fits'
    dateobs = product.obs_beg.strftime("%Y%m%dT%H%M%S")
    dateend=product.obs_end.strftime("%Y%m%dT%H%M%S")
    if run_type=='file': #to be consistent with the old convention
        run_type=''
    else:
        run_type='_'+run_type
    
    return f'solo_L1A_stix-{product_type.replace("_", "-")}{run_type}_{dateobs}-{dateend}_{unique_id:06d}_V{version:02d}.fits'
            #hk-mini or hk-maxi


def write_fits(basepath,unique_id, prod, product_type, overwrite=True, version=1, run_type='file'): 
    filename = generate_filename(prod, unique_id, product_type, version=version, run_type=run_type)
    primary_header = generate_primary_header(filename, prod.scet_coarse, prod.scet_fine,
                                                 prod.obs_beg, prod.obs_avg, prod.obs_end, version)
    hdul = prod.to_hdul()
    hdul[0].header.update(primary_header)
    hdul[0].header.update({'HISTORY': 'Processed by STIX LLDP VM'})

    full_path=basepath/filename
    if full_path.is_file():
        print("Removing existing fits:", str(full_path))
        full_path.unlink()
        

    hdul.writeto(full_path, checksum=True, overwrite=overwrite)
    return {'data_start_unix': prod.obs_beg.timestamp(), 'data_end_unix':prod.obs_end.timestamp(), '_id':unique_id,
            'filename': filename}

