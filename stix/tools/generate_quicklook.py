import pymongo
import sys
from stix.core import datatypes as sdt
from stix.spice import datetime

class StixQuickLookReportAnalyzer(object):
    """
    capture quicklook reports and fill packet information into a MongoDB collection

    """
    def __init__(self, collection):
        self.db_collection = collection
        self.report = None
        try:
            self.current_report_id = self.db_collection.find().sort(
                '_id', -1).limit(1)[0]['_id'] + 1
        except IndexError:
            self.current_report_id = 0
        self.last_run=0
    

    def capture(self, pkt):
        run_id=pkt['run_id']
        if self.last_run!=run_id:
            print(run_id)
        self.last_run=run_id
        packet_id=pkt['_id']
        if not self.db_collection:
            return
        packet = sdt.Packet(pkt)

        if not packet['parameters']:
            return


        if packet.SPID not in [54118, 54119, 54121,54120]:
            return

        start_coarse_time=0
        start_fine_time=0
        integrations=0
        detector_mask=0
        pixel_mask=0
        start_unix_time=0
        duration=0
        start_coarse_time = packet[1].raw
        start_fine_time = packet[2].raw
        integrations = packet[3].raw
        points=0

        if packet.SPID == 54118:
            #QL LC
            #detector_mask = packet[4].raw
            #pixel_mask = packet[6].raw
            points = packet.get('NIX00270/NIX00271')[0][0]
        elif packet.SPID==54119:
            points=packet[15].raw
        elif packet.SPID==54121:
            points=packet[13].raw
        elif packet.SPID==54120:
            points=packet[14].raw


        duration = points * 0.1 * (integrations + 1)
        start_unix_time = datetime.scet2unix(start_coarse_time,
                                                      start_fine_time)



        report={
            '_id': self.current_report_id,
            'run_id': run_id,
            'SPID': packet.SPID,
            'packet_id': packet_id,
            'start_unix_time': start_unix_time,
            'packet_header_time': packet['unix_time'],
            'integrations': integrations,
            'duration': duration,
            'stop_unix_time': start_unix_time + duration,
            #'detector_mask': detector_mask,
            #'pixel_mask': pixel_mask
        }
        self.db_collection.insert_one(report)
        self.current_report_id += 1

connect = pymongo.MongoClient('localhost', 27017)

stix=connect['stix']

ql=StixQuickLookReportAnalyzer(stix['quick_look'])
db=stix['packets']
cur=db.find({'header.SPID':{'$in':[54118,54119,54120,54121]}, 'run_id':{'$gt':113}}).sort('_id',1)

for pkt in cur:
    ql.capture(pkt)



