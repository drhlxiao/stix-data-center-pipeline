#!/usr/bin/python
#Fetch fermi triggers
#Hualin Xiao Dec. 16 2016

import sys

from log_file import get_logfile
import re
import os
import math
import time
from dateutil import parser
from urllib.request import urlopen
from urllib.error import URLError
from argparse import ArgumentParser
from bs4 import BeautifulSoup
from datetime import datetime

from astropy import units as u
from astropy.coordinates import SkyCoord
from  grb_db import *
import logging



def convert_ra_dec(rah, decd):
    try:
        rah=rah.strip()
        decd=decd.strip()
        coor=rah+" "+decd
        c = SkyCoord(coor, unit=(u.hourangle, u.deg))
        coordeg=[c.ra.degree, c.dec.degree]
        return coordeg
    except:
        return []



def insert_database(trigid, gdatetime, ra, dec,url):
    
    msg=""
    if not grb_exists('INTEGRAL',trigid, gdatetime):
        logging.info( "inserting event:")
        logging.info( "               : "+trigid+" "+gdatetime)
        ins="insert into burst_events(name, start_utc, stop_utc, burst_utc, right_ascension, declination, url, comments, author) values('INTEGRAL %s', '%s', '%s', '%s', %f,%f,'%s','Also detected by INTEGRAL','robot')"%(trigid, gdatetime,gdatetime,gdatetime,float(ra),float(dec),url)
        print(ins)
        logging.info( ins)
        excute_sql(ins)
        msg="Integral %s    %s   %f   %f"%(trigid, gdatetime, float(ra),float(dec))
    else:
        logging.info( "Already existed Integral trigger:")
        logging.info( "                   :%s \t %s \t %s \t %s"%(trigid, gdatetime,ra, dec))

    
    return msg

def extract_LC_url(item):
    iurl=""
    for link in item.find_all('a'):
        iurl=link.get('href')
    return  iurl 
def is_grb_name(name):
    name2=name.strip()
    grbname= re.findall('\d{6}[A-Z]', name2)
    if grbname:
        return True
    else:
        return False




def parse_rows(rows):
    """ Get data from rows """
    msglist= []
    nn=0
    for row in rows:
        table_data = row.find_all('th')
        if len(table_data)>1:
            grb=table_data[0].get_text()
            if is_grb_name(grb):
                try:
                    imgurl=extract_LC_url(table_data[11])
                    rah=table_data[1].get_text()
                    decd=table_data[2].get_text()
                    angles=convert_ra_dec(rah,decd)
                    if len(angles)==0:
                        continue
                    ra=angles[0]
                    dec=angles[1]
                    gdate=grb[0:7].strip()
                    gtime=table_data[5].get_text()
                    gdatetime=datetime.strptime(gdate,'%y%m%d').strftime("%Y-%m-%d")+" "+gtime
                    integral_grbname="GRB "+grb 
                    msg=insert_database(integral_grbname, gdatetime, ra, dec,imgurl)

                    nn=nn+1
                except:
                       continue
            else:
                logging.info( 'not a valid grb name')



def main():
    # Get arguments
    # Make soup
    url='http://www.isdc.unige.ch/integral/ibas/cgi-bin/ibas_isgri_web.cgi'

    logging.info( 'fetching web page...')
    try:
        resp = urlopen(url)
    except URLError as e:
        logging.info( 'An error occured fetching %s \n %s' % (url, e.reason)   )
        return 1
    content=f"<html><body>{resp.read()}</body></html>"
    soup = BeautifulSoup(content,"html.parser")
    # Get table
    logging.info( 'finding table...')
    try:
        table = soup.find('table',{"class":"isdctable"})
    except AttributeError as e:
        logging.error( 'No tables found, exiting')
        return 1

    # Get rows
    logging.info( 'finding rows...')
    try:
        rows = table.find_all('tr')
    except AttributeError as e:
        logging.error( 'No table rows found, exiting')
        return 1

    logging.info( 'parsing table...')
    table_data = parse_rows(rows)


if __name__ == '__main__':

    
    logfile=get_logfile()
    logging.basicConfig(level=logging.DEBUG,
                     format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=logfile,
                    filemode='a')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    logging.info('----------------------------------------------')


    status = main()
    sys.exit(status)
