import grb_db

from dateutil import parser as dtparser
import os
import sys
sys.path.append('.')
sys.path.append('../')
sys.path.append('../../')
from scipy import signal
import numpy as np
import math
import matplotlib
from matplotlib import pyplot as plt
from stix.core import stix_datatypes as sdt
from stix.core import mongo_db as db
from stix.spice import stix_datetime
from stix.core import stix_logger
from stix.spice import solo
from stix.flare_pipeline import stix_lightcurves as stl
from stix.analysis import ql_analyzer as qla

def main():
    sql='select start_utc, name  from  burst_events'
    rows=grb_db.fetchall(sql)
    urls=[]
    for row in rows:
        print(row[0])
        start=sdt.utc2unix(row[0])-900
        end=start+1800
        su=stix_datetime.unix2utc(start)
        eu=stix_datetime.unix2utc(end)
        print(row)
        try:
            stl.plot_stix_lc('LC_GRBs', 0, f'{row[0]}_{row[1]}', su, eu, True, event_type=' ')
        except Exception as e:
            print(str(e))



if __name__=='__main__':
    main()
