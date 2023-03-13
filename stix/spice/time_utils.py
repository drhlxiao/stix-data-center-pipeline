#!/usr/bin/python3

import re
import glob
from datetime import datetime
from dateutil import parser as dtparser
from astropy.time import Time
from stix.spice import spice_manager as spm
import spiceypy
from datetime import datetime
import pandas as pd
from astropy.time import Time

from stix.core import logger

def anytime(dt, fm='iso'):
    if isinstance(dt, Time):
        dt=dt.to_datetime()
    t=pd.to_datetime(dt,utc=True)
    if fm=='iso':
        return t.strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    elif fm=='unix':
        return t.timestamp()
    elif fm=='datetime':
        return t.to_pydatetime()
    return t



def utc2unix(dt):
    t=pd.to_datetime(dt,utc=True)
    return t.timestamp()
    

def utc2datetime(t):
    return pd.to_datetime(t, utc=True).to_pydatetime()

def datetime2unix(t):
   t=pd.to_datetime(t, utc=True)
   return t.timestamp()




def unix2datetime(unix_timestamp):
    return datetime.utcfromtimestamp(unix_timestamp)

def format_datetime(dt):
    if isinstance(dt, datetime):
        return dt.isoformat(timespec='milliseconds')
    elif isinstance(dt, (int, float)):
        return datetime.utcfromtimestamp(dt).isoformat(timespec='milliseconds')
    elif isinstance(dt, str):
        try:
            return format_datetime(float(dt))
        except ValueError:
            return dt
    else:
        return '1970-01-01T00:00:00.000Z'

def now():
    return datetime.utcnow()

def get_now(dtype='unix'):
    utc_iso = datetime.utcnow().isoformat() + 'Z'
    if dtype == 'unix':
        return dtparser.parse(utc_iso).timestamp()
    return utc_iso


def scet2utc(coarse, fine=0):
    try:
        coarse_int = coarse
        if isinstance(coarse, float):
            coarse_int = int(coarse)
        fine_int = int((coarse - coarse_int) * 65536)+fine

        return spm.spice.scet2utc(coarse_int, fine_int)
    except spiceypy.utils.support_types.SpiceyError:
        return ''


def utc2scet(utc):
    return spm.spice.utc2obt(utc)



def utc2isoformat(utc):
    if not utc.endswith('Z'):
        utc += 'Z'
    try:
        return dtparser.parse(utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    except:
        return None



def unix2datetime(timestamp):
    dt = None
    if isinstance(timestamp, float):
        dt = datetime.utcfromtimestamp(timestamp)
    elif isinstance(timestamp, str):
        try:
            ts = float(timestamp)
            dt = datetime.utcfromtimestamp(ts)
        except ValueError:
            dt = dtparser.parse(timestamp)
    elif isinstance(timestamp, datetime.datetime):
        dt = timestamp
    return dt



def scet2unix(coarse, fine=0):
    try:
        utc = scet2utc(coarse, fine)
        return utc2unix(utc)
    except spiceypy.utils.support_types.SpiceyError:
        return 0


def unix2utc(ts):
    return datetime.utcfromtimestamp(ts).isoformat(timespec='milliseconds')


def unix2scet(unix):
    utc = unix2utc(int(unix))
    return utc2scet(utc)

def utc2filename(utc):
    dt=utc2datetime(utc)
    return dt.strftime("%Y%m%dT%H%M%S")

def scet2datetime(coarse, fine=0):
    unixtimestamp = scet2unix(coarse, fine)
    return datetime.utcfromtimestamp(unixtimestamp)


def unix2datetime(unix_timestamp):
    return datetime.utcfromtimestamp(unix_timestamp)

def is_unix_time_valid(unix):
    unix=float(unix)
    return unix > 1577833200 and unix< 2524604400
    #2020-01-01 to 2050-01-01
def is_scet_valid(scet):
    j2000=946724400
    scet=float(scet)
    return is_unix_time_valid(j2000+scet)


'''
    code taken from Shane's datetime script
'''






def scet_to_utc(scet):
    return spm.spice.obt2utc(scet)


def scet_to_datetime(scet):
    utc_iso = scet_to_utc(scet)
    return dtparser.isoparse(utc_iso)


def utc_to_scet(utc):
    return spm.spice.utc2scet(utc)


def datetime_to_scet(dt):
    if isinstance(dt, Time):
        dt = dt.to_datetime()
    utc_iso = dt.isoformat(timespec='microseconds')
    scet = utc_to_scet(utc_iso)[2:]
    return scet


if __name__ == '__main__':
    print("UTC at T0:")
    print(scet2utc(0))
