#!/usr/bin/python
"""
Image viewer, create image from fits files which are created by idl imaging software 

April 27, 2022
"""
import os
import matplotlib
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib import pyplot as plt
from dateutil.parser import parse as dtparse

from sunpy import map
from sunpy.map import make_fitswcs_header

from sunpy.coordinates.frames import HeliocentricEarthEcliptic,HeliographicStonyhurst

from stix.flare_pipeline import lightcurves
from stix.core import mongo_db as db

mdb = db.MongoDB()
flare_images_db= mdb.get_collection('flare_images')

SMALL_SIZE = 8
matplotlib.rc('font', size=SMALL_SIZE)
matplotlib.rc('axes', titlesize=SMALL_SIZE)
matplotlib.rcParams['axes.titlepad']=35

CMAP='RdPu' #color map

def create_images_in_queue():
    cursor=flare_images_db.find({'num_idl_calls':0, 'idl_status':''}).sort('_id',-1)
    create_images_for_bsd_docs(cursor)

def create_images_for_science_data(bsd_id):
    cursor=flare_images_db.find({'bsd_id':bsd_id})
    create_images_for_bsd_docs(cursor)

def create_qk_figures_for_all():
    docs=flare_images_db.find({'idl_status':True, 'figs.0':{'$exists':False}}).sort('_id',-1)
    for doc in docs:
        create_qk_figures(doc, update_db=True)
def process_one_latest():
    cursor=flare_images_db.find({'idl_status':'', 'figs.0':{'$exists':False}}).sort('_id',-1).limit(1)
    for doc in cursor:
        logger.info("Processing doc: {doc['_id']}")
        process_one(doc)


def create_images_for_bsd_docs(cursor):
    for doc in cursor:
        process_one(doc)

def create_qk_figures(doc, update_db=False, output_folder=None):
    """
    write images to the same folder if output_folder is None
    """
    try:
        solo_hee=np.array(doc['aux']['solo_hee'][0])*u.km
        fnames=[doc['fits'][0], doc['fits'][1]]
        start_utc, end_utc=doc['utc_range']
        energy_range=doc['energy_range']
        logger.info(f'Creating images for flare image #{doc["_id"]}')
        figs=imv.images_to_graph(fnames[0], fnames[1], solo_hee, doc['bsd_id'], 
                start_utc, end_utc, energy_range, output_folder, doc['unique_id'], doc['background']['unique_id'] )
        logger.info(str(figs))
    except Exception as e:
        raise
        logger.error(e)
        return None
    if update_db:
        updates={'$set':{'figs':figs, 'processing_date':datetime.now()}}

        flare_images_db.update_one({'_id':doc['_id']}, updates)
    return figs

def create_figures_ids_between(start_id, end_id):
    for i in range(start_id, end_id):
        doc=flare_images_db.find_one({'_id':i})
        if not doc:
            logger.warning(f"Failed to create figures for DocID:{i}")
            continue
        logger.info(f"Creating images for BSD#{doc['bsd_id']}, DocID:{i}")
        create_qk_figures(doc, update_db=True, output_folder=None)
def process_one(doc):
    args=doc.get('idl_args',None)
    if args is not None:
        logger.info(f"Processing {doc['_id']}")
        logger.info(f": prameters {str(args)}")
        success = call_idl(args[0], args[1], args[2], args[3])
        logger.info(f"End of processing {doc['_id']}, status: {success}")
        updates={'num_idl_calls':doc['num_idl_calls']+1, 'idl_status':success}
        if success:
            figs=create_qk_figures(doc, update_db=None, output_folder=None)
            if figs:
                updates['figs']=figs
                updates['processing_date']=datetime.now()
        updates={'$set':updates}
        flare_images_db.update_one({'_id':doc['_id']}, updates)

            

def create_STIX_map(fits_filename,  solo_hee,   map_name=''): 

    """
    Converts the input FITS to a standard Sunpy Solar map, by
    returning the map containing the STIX source.

    Parameters
    ----
    fits_filename: str
        FITS filename 
    solo_hee:  Quantity(u.km)
        solar orbiter coordinates in the HeliocentricEarthEcliptic reference frame
    map_name: str, optional
        map name, 
        
    Returns
    ---
    stix_map: sunpy.map
        STIX map
    indices:  np.array
        pixel indices at maximum
    header:  dict
        fits file primary HDUlist header

    """
    hdul=fits.open(fits_filename)
    header=hdul[0].header
    obs_time=header['DATE_OBS']
    
    # Convert the string date and time to datetime.datetime
    datetime_map = dtparse(obs_time)
    # Obtain the HEE coordinates of Solar Orbiter

    solo_hee = HeliocentricEarthEcliptic(solo_hee,
                                              obstime=Time(datetime_map).isot,
                                              representation_type='cartesian')
    image_center=[0,0];
    image_center[0]=header['NAXIS1']/2.+0.5 - header['CRPIX1'] #
    image_center[1]=header['NAXIS2']/2.+0.5 - header['CRPIX2']
    print(image_center)
    # Set the coordinates of the reference pixel of the STIX map
    stix_ref_coord = SkyCoord(image_center[0]*u.arcsec, 
                              image_center[1]*u.arcsec,
                              obstime=datetime_map,
                              observer=solo_hee.transform_to(HeliographicStonyhurst(obstime=datetime_map)),
                              frame='helioprojective')
    
    # Get the distance of Solar Orbiter from the Sun (center)

    # Create a FITS header containing the World Coordinate System (WCS) information
    scale=np.array([header['CDELT1'], header['CDELT2']])
    

    out_header = make_fitswcs_header(hdul[0].data,
                                     coordinate=stix_ref_coord, 
                                     #reference_pixel=np.array([header['CRPIX1'], header['CRPIX2']])*u.pixel,
                                     
                                     scale=scale*u.arcsec/u.pixel,
                                     rotation_angle =-header['CROTA1'] * u.deg,
                                     instrument=f"STIX {map_name} map",
                                     observatory="Solar Orbiter")
    max_value=hdul[0].data.max()
    max_indices = np.where(hdul[0].data ==max_value)



    stix_map = map.Map(hdul[0].data, out_header)
    return stix_map, max_indices, header
def images_to_graph(fw_image_fname, bp_image_fname, solo_hee, bsd_id,  start_utc, end_utc, energy_range,  output_folder=None, bsd_uid=None, bkg_uid=None, processing_id=None):
    """
    Plot STIX image
    Parameters
    ----------
    bp_image_fname: str
        filename of fits file containing the Back-projection image
    fw_image_fname: str
        filename of fits file containing the forward-fit image    
     solo_hee:  Quantity(u.km)
        solar orbiter coordinates in the HeliocentricEarthEcliptic reference frame
    output_fname_prefix: str
        output filename
    Returns
    -------
    ret:  dict
         json string to reproduce figures
    bmap: sunpy.map
        map of back-projection image
    fwmap: sunpy.map
        map of forward fit image

    """
    print("args:", bp_image_fname, fw_image_fname, solo_hee, bsd_id, start_utc, end_utc, energy_range,  output_folder, "end")
    if not os.path.isfile(bp_image_fname) or not os.path.isfile(fw_image_fname):
        print(bp_image_fname, fw_image_fname)
        return {'error': 'Image files not found'}



    map_name=''
    ret={}
    bp_map, bp_maxidx,bp_header=create_STIX_map(bp_image_fname,   solo_hee,
                                            map_name)
    
    fw_map, fw_maxidx, fw_header=create_STIX_map(fw_image_fname,    
                                                solo_hee,map_name)

    title = f"{bp_header['ORIGIN']} \n {bp_header['DATE_OBS']} "

    fig = plt.figure(figsize=(9,7), dpi=200, facecolor='white')

    ax_lc= fig.add_subplot(221)
    #light curves
    lightcurves.plot_QL_lc_for_bsd(bsd_id, fill_between_times=[start_utc, end_utc], ax=ax_lc)

    #back projection
    ax_bp = fig.add_subplot(222, projection=bp_map)
    bp_map.plot( vmin=bp_map.max()*0.1, 
                    cmap=CMAP,
                    axes=ax_bp,  title="Back-projection",)

    bp_map.draw_grid(color='k', ls='--', grid_spacing=10*u.deg)
    bp_map.draw_limb(axes=ax_bp, color='k',alpha=0.5)



    ax_bp_clicp = fig.add_subplot(223, projection=bp_map)
    bp_map.plot( vmin=bp_map.max()*0.1, 
                    cmap=CMAP,
                    axes=ax_bp_clicp, title='Back-projection clipped')
    bp_map.draw_limb(axes=ax_bp_clicp, color='k',alpha=0.5)

    ymin,ymax=ax_bp.get_ylim()
    xmin,xmax=ax_bp.get_xlim()

    print("axis ranges", xmin, xmax, ymin, ymax)
    hw,hh=50, 50

    xlim=(max(xmin,bp_maxidx[1][0]-hw), min(xmax, bp_maxidx[1][0]+hw))
    ylim=(max(ymin,bp_maxidx[0][0]-hh), min(ymax, bp_maxidx[0][0]+hh))

    print("lims", xlim, ylim)
    try:
        ax_bp_clicp.set_xlim(xlim)
        ax_bp_clicp.set_ylim(ylim)
    except Exception as e:
        logger.warning(f"Failed to clip image: {e}")
        #astropy may raise an error when writing the image to a file, we need to modify the astropy
        #find the line change it to the line as below in visualization/wcsaxes/transforms.py 
        #if self.same_frames or not np.any(input_coords):
        #    return input_coords
        pass

    bp_map.draw_grid(color='k', ls='--', grid_spacing=10*u.deg)
    #forward fit

    ax_fw = fig.add_subplot(224, projection=fw_map)
    fw_map.plot(vmin=0., 
                            cmap=CMAP,
                            axes=ax_fw, title='Forward-fit')
    fw_map.draw_limb(axes=ax_fw, color='k', alpha=0.5)
    fw_map.draw_grid(color='k', ls='--', grid_spacing=10*u.deg)


    if output_folder is None:
        output_fname_prefix=bp_image_fname.rsplit('.', 1)[0] 
    else:
        basename=os.path.basename(bp_image_fname)
        prefix=os.path.join(output_folder, basename)
        output_fname_prefix=prefix.rsplit('.', 1)[0] 

    output_fname=f'{output_fname_prefix}.svg' 
    process_id_str = f'(#{processing_id})' if processing_id is not None else ''
    plt.suptitle(f'STIX light curves and images {process_id_str} \n IMG EXP TIME {start_utc} –  {end_utc}\n E: {energy_range[0]} – {energy_range[1]} keV  SIG/BKG: {bsd_uid}/{bkg_uid} ')
    #plt.tight_layout()
    plt.subplots_adjust(
            top=.82,
                    wspace=0.1, 
                    hspace=0.6)
    plt.savefig(output_fname, format='svg')
    print(output_fname)

    return [output_fname]



def test_image_viewer():
    f1='/data/quicklook/flare_images/sci_9111_uid_2202155513_0_0_7_2022-02-15T17:21:49.918_bp.fits'
    f2='/data/quicklook/flare_images/sci_9111_uid_2202155513_0_0_7_2022-02-15T17:21:49.918_fwfit.fits'
    hee=np.array([ 1.03180004e+08,  -3.22106960e+07,  5.13784704e+06])*u.km
    images_to_graph(f2, f1, hee, 9111, '2022-02-15T17:21:49.918', '2022-02-15T17:21:51.918', [4,10], '.', 0, 0 )
