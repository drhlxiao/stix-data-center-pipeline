
import os
import sys
import numpy as np
from datetime import datetime
from stix.spice import stix_datetime
from stix.core import stix_datatypes as sdt
from stix.core import mongo_db as mgdb
from stix.core import stix_metadata as md
mdb = mgdb.MongoDB()

QLLC_SPID = 54118
QLBKG_SPID= 54120
db_pkt=mdb.get_collection('packets')
db=mdb.get_db()
id_start=812

packets=db_pkt.find({'header.SPID':54122,'run_id':{'$gt':id_start}}).sort('_id',1)

qla= md.StixQuickLookReportAnalyzer(db)
last_run=0

for pkt in packets:
    run_id=pkt['run_id']
    packet_id=pkt['_id']
    if last_run!=run_id:
        print(run_id)
    last_run=run_id
    if pkt['header']['SPID']==54122:
        qla.write_ql_flare_loc_to_db(pkt, run_id,packet_id)
    #if pkt['header']['SPID']==54120:
    #    qla.write_ql_spec_to_db(pkt, run_id,packet_id)
    #else:
    #    qla.write_qllc_to_db(pkt, run_id,packet_id)

