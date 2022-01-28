import sys
import os
import pprint
dirpath='../parser'
idb_fname=dirpath+'/idb/idb.sqlite'
from core import idb
idb = idb.idb(idb_fname)
pprint.pprint(idb.get_fixed_packet_structure(54102))
def get_fixed_telemetry_length(spid):
    rows=idb.get_fixed_packet_structure(int(spid))
    if rows:
        last_row=rows[-1]
        return last_row['PLF_OFFBY']+last_row['PCF_WIDTH']/8
    else:
        return None

def get_telemetry_length():
    sql='select PID_SPID, PID_APID, PID_TYPE,PID_STYPE, PID_DESCR, PID_TPSD from PID  order by PID_TPSD asc'
    rows = idb.execute(sql, None, 'list')
    for res in rows:
        spid, apid, ser_type,ser_stype, desc, tpsd=res
        length='variable'
        if tpsd==-1:
            length=get_fixed_telemetry_length(spid)
        print('{}; {};{};{};{};{}'.format(spid, apid, ser_type,ser_stype, length, desc))

get_telemetry_length()
