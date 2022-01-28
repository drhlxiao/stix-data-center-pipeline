import sys
import os
from collections import defaultdict
from datetime import datetime
from itertools import chain
from pathlib import Path

from stix.core import datatypes as sdt
from stix.fits.io.processors import FitsL1Processor
from stix.fits.io import hk_fits_writer as hkw
from stix.spice.datetime import datetime_to_scet
from stix.fits.products.housekeeping import MiniReport, MaxiReport
from stix.fits.products.quicklook import LightCurve, Background, Spectra, Variance, \
    FlareFlagAndLocation, CalibrationSpectra, TMManagementAndFlareList
from stix.fits.products.science import XrayL0, Aspect, XrayL1, XrayL2, XrayL3, Spectrogram
#from stix.utils.logger import get_logger


from stix.core import mongo_db, logger
logger = logger.get_logger()


db= mongo_db.MongoDB()
#logger = get_logger(__name__)

FITS_PATH='/data/fits/'

DATA_LEVEL='L1A'



SPID_MAP = {
    # Bulk Science
    #54110: 'xray_l0_auto',
    #54111: 'xray_l1_auto',
    #54112: 'xray_l2_auto',
    #54113: 'xray_l3_auto',
    #54142: 'spectrogram_auto',
    #54114: 'xray_l0_user',
    #54115: 'xray_l1_user',
    #54116: 'xray_l2_user',
    #54117: 'xray_l3_user',
    #54143: 'spectrogram_user',
    54125: 'aspect',
    # Quick look
    #54118: 'ql_light_curves',
    #54119: 'ql_background',
    #54120: 'ql_spectrogram',
    #54121: 'ql_variance',
    #54122: 'flareflag_location',
    #54123: 'tm_status_and_flare_list',
    #54124: 'calibration_spectrum',
    # House keeping
    #54101: 'hk_mini',
    #54102: 'hk_maxi'
}
SEG_FLAG_MAP={0: 'continuation packet',1: 'first packet',2: 'last_packet',3:'stand-alone packet'}



def process_packets(file_id, packet_lists, spid, product, report_status,  base_path_name=FITS_PATH, overwrite=True, version=1):
    """
    Process a sequence containing one or more packets for a given product.

    Parameters
    ----------
    packet_lists : list
        Packets
    spid : int
        SPID
    product : basestring
        Product name
    basepath : pathlib.Path
        Path
    pacekt_type:
        complete packets or incomplete packets
    overwrite : bool (optional)
        False (default) will raise error if fits file exits, True overwrite existing file
    """

    """
    Process a sequence containing one or more packets for a given product.

    Parameters
    ----------
    packet_lists : list
        List of packet sequences
    spid : int
        SPID
    product : basestring
        Product name
    basepath : pathlib.Path
        Path
    overwrite : bool (optional)
        False (default) will raise error if fits file exits, True overwrite existing file
    """
    # For HK merge all stand alone packets in request
    if spid in [54101, 54102]:
        packet_lists = [list(chain.from_iterable(packet_lists))]
    #print('BASE',base_path_name)

    for packets in packet_lists:
        if not packets:
            continue
        parsed_packets = sdt.Packet.merge(packets, spid, value_type='raw')
        e_parsed_packets = sdt.Packet.merge(packets, spid, value_type='eng')


        try:
            if product == 'hk_mini':
                prod = MiniReport(parsed_packets)
                product_type = 'housekeeping'
            elif product == 'hk_maxi':
                prod = MaxiReport(parsed_packets)
                product_type = 'housekeeping'
            elif product == 'ql_light_curves':
                #prod = LightCurve(parsed_packets, e_parsed_packets)
                prod = LightCurve.from_packets(parsed_packets, e_parsed_packets)
                product_type = 'quicklook'
            elif product == 'ql_background':
                prod = Background.from_packets(parsed_packets, e_parsed_packets)
                #prod = Background(parsed_packets, e_parsed_packets)
                product_type = 'quicklook'
            elif product == 'ql_spectrogram':
                prod = Spectra.from_packets(parsed_packets, e_parsed_packets)
                #prod = Spectra(parsed_packets, e_parsed_packets)
                product_type = 'quicklook'
            elif product == 'ql_variance':
                prod = Variance.from_packets(parsed_packets, e_parsed_packets)
                product_type = 'quicklook'
            elif product == 'flareflag_location':
                prod = FlareFlagAndLocation.from_packets(parsed_packets)
                product_type = 'quicklook'
            elif product == 'calibration_spectrum':
                prod = CalibrationSpectra.from_packets(parsed_packets, e_parsed_packets)
                product_type = 'quicklook'
            elif product == 'tm_status_and_flare_list':
                prod = TMManagementAndFlareList.from_packets(parsed_packets)
                product_type = 'quicklook'
            elif product == 'xray_l0_user':
                prod = XrayL0.from_packets(parsed_packets, e_parsed_packets)
                product_type = 'science'
            elif product == 'xray_l1_user':
                prod = XrayL1.from_packets(parsed_packets, e_parsed_packets)
                product_type = 'science'
            elif product == 'xray_l2_user':
                prod = XrayL2.from_packets(parsed_packets, e_parsed_packets)
                product_type = 'science'
            elif product == 'xray_l3_user':
                prod = XrayL3.from_packets(parsed_packets, e_parsed_packets)
                product_type = 'science'
            elif product == 'spectrogram_user':
                prod = Spectrogram.from_packets(parsed_packets, e_parsed_packets)
                product_type = 'science'
            elif product == 'aspect':
                prod = Aspect.from_packets(parsed_packets, e_parsed_packets)
                product_type = 'science'
            else:
                logger.warning(f'Not implemented {product}, SPID {spid}.')
                continue



            base_path=Path(base_path_name)
            base_path.mkdir(parents=True, exist_ok=True)

            unique_id=db.get_next_fits_id()
            #print('Unique id:',unique_id)
            complete=report_status=='complete'
            meta_entries=[]
            if product_type=='housekeeping':
                meta=hkw.write_fits(base_path,unique_id,  prod, product, overwrite, version) 
                meta_entries=[meta]
            else:
                fits_processor = FitsL1Processor(base_path, unique_id, version)
                fits_processor.write_fits(prod)
                meta_entries=fits_processor.get_meta_data()
            for meta in meta_entries :
                file_size=0
                try: 
                    abs_path=base_path/meta['filename']
                    file_size=abs_path.stat().st_size
                except:
                    raise
                    pass

                doc={
                    #'_id':unique_id,
                    'packet_id_start': packets[0]['_id'],
                    'packet_id_end': packets[-1]['_id'],
                    'packet_spid':spid,
                    'num_packets': len(packets),
                    'file_id':int(file_id), 
                    'product_type':product, 
                    'product_group':product_type,
                    #'data_start_unix':meta['data_start_unix'],
                    #'data_end_unix':meta['data_end_unix'],
                    #'filename': meta['filename'],
                    'complete':complete,
                    'version': version,
                    'level':DATA_LEVEL,
                    'creation_time':datetime.utcnow(),
                    'path':base_path_name,
                    'file_size':file_size
                    }
                doc.update(meta)
                #print(doc)
                db.write_fits_index_info(doc)
                logger.info(f'created  fits file:  {meta["filename"]}')

        except Exception as e:
            raise
            logger.error(str(e))

def purge_fits_for_raw_file(file_id):
    print(f'Removing existing fits files for {file_id}')
    fits_collection=db.get_collection('fits')
    if fits_collection:
        cursor=fits_collection.find({'file_id':int(file_id),'product_type':'aspect'})
        for cur in cursor:
            try:
                fits_filename=os.path.join(cur['path'],cur['filename'])
                logger.info(f'Removing file: {fits_filename}')
                os.unlink(fits_filename)
            except Exception as e:
                raise
                logger.warning(f'Failed to remove fits file:{fits_filename} due to: {str(e)}')
        cursor = fits_collection.delete_many({'file_id': int(file_id), 'product_type':'aspect'})






def create_fits(file_id, output_path, overwrite=True,  version=1):
    if overwrite:
        purge_fits_for_raw_file(file_id)
    spid_packets=db.get_file_spids(file_id)
    #print(spid_packets)
    if not spid_packets:
        logger.warning(f'File {file_id} has no packet!')
        return
    for spid in spid_packets:
        spid=int(spid)
        if spid not in SPID_MAP.keys():
            logger.warning(f'Not supported spid : {spid}')
            continue
        logger.info(f'Requesting packets of file {file_id} from MongoDB')
        sort_field='header.unix_time'
        #if spid==54118:
        cursor= db.select_packets_by_run(file_id, [spid],sort_field=sort_field)
        if cursor:
            logger.info(f'packets have been selected for SPID {spid}')

        product = SPID_MAP[spid]
        report_status='complete'
        if spid in [54101, 54102]:
            packets=[list(cursor)]
            process_packets(file_id, packets, spid, product, report_status, output_path, overwrite, version)
            continue
        packets=[]
        for i, pkt in enumerate(cursor):
            seg_flag = int(pkt['header']['seg_flag'])
            #print(f'Packet: {i} - { SEG_FLAG_MAP[seg_flag]}')
            
            if seg_flag == 3: #'stand-alone packet':
                packets=[pkt]
            elif seg_flag == 1:# 'first packet':
                packets= [pkt]
            elif seg_flag == 0: #'continuation packet':
                if packets:
                    packets.append(pkt)
                else:
                    packets= [pkt]
            elif seg_flag ==2: # 'last packet':
                if packets:
                    packets.append(pkt)
                else:
                    packets= [pkt]
            if seg_flag in [3,2]:
                report_status='complete'
                #print('complete packets')
                process_packets(file_id, [packets], spid, product, report_status, output_path, overwrite, version)
                packets=[] #clean container after processing packets
        if packets:
            #incomplete packets, QL packets are always incomplete 
            report_status='incomplete'
            #print('incomplete packets:')
            #print([(pkt['_id'], pkt['header']['segmentation']) for pkt in packets])
            process_packets(file_id, [packets], spid, product, report_status, output_path, overwrite, version)




if  __name__ == '__main__':
    if len(sys.argv) == 1:
        print('''STIX L1 FITS writer
                   It reads packets from mongodb and write them into fits
                Usage:
                packets_to_fits   <file_id>
                packets_to_fits   <file_start_id> <end_id>
                ''')
    elif len(sys.argv)==2:
        create_fits(sys.argv[1], FITS_PATH, overwrite=True, version=1)
    else:
        start=int(sys.argv[1])
        end=int(sys.argv[2])
        for i in range(start,end+1):
            purge_fits_for_raw_file(i)
        for i in range(start,end+1):
            create_fits(i, FITS_PATH, overwrite=True, version=1)

