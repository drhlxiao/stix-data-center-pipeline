#!/usr/bin/python

#Fetch knous triggers
#Hualin Xiao Dec. 16 2016

import sys

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
from  grb_db import *
import logging

from log_file import get_logfile





def insert_database(trigid, gdatetime, ra, dec,url):
    
    msg=""
    if not grb_exists('KONUS',trigid, gdatetime):
        logging.info( "inserting...")
        ins="insert into burst_events(name, start_utc, stop_utc, burst_utc, right_ascension, declination, url, comments, author) values('KONUS %s', '%s', '%s','%s', %f,%f,'%s','detected by Konus(%s)','robot')"%(trigid, gdatetime,gdatetime,gdatetime,float(ra),float(dec),url,trigid)
        logging.info( ins)
        excute_sql(ins)
        msg="GBM %s    %s   %f   %f"%(trigid, gdatetime, float(ra),float(dec))
    else:
        logging.info( "Already existed fermi trigger:")
        logging.info( "                   :%s \t %s \t %s \t %s"%(trigid, gdatetime,ra, dec))

    
    return msg


def parse_rows(rows):
    """ Get data from rows """
    msglist= []
    nn=0
    for row in rows:
        table_data = row.find_all('td')
        if len(table_data)>10:
            dd=str(table_data[0].get_text())
            if nn>10:
                break  #don't process more than this number per day

            if dd:
                try:
                    dt=table_data[2].get_text()
                    lc=table_data[7]
                    gtime=re.findall(r"\d{2}:\d{2}:\d{2}",dt)[0]
                    

                    gdatetime=datetime.strptime(str(dd).strip(),'%Y%m%d').strftime("%Y-%m-%d")+" "+gtime
                    logging.info( gdatetime)
                    link=lc.find_all('a')[0]
                    imgurl="http://gcn.gsfc.nasa.gov/"+link.get('href')
                    msg=insert_database(dd,gdatetime,-1,-1,imgurl)
                    nn=nn+1
                except:
                    continue




def main():
    # Get arguments
    # Make soup
    url='https://gcn.gsfc.nasa.gov/konus_grbs.html'
    try:
        logging.info( "opening %s"%url)
        resp = urlopen(url)
    except URLError as e:
        logging.info( 'An error occured fetching %s \n %s' % (url, e.reason))
        return 1
    soup = BeautifulSoup(resp.read(),"html.parser")
    # Get table
    try:
        table = soup.find('table')
    except AttributeError as e:
        logging.info( 'No tables found, exiting')
        return 1

    # Get rows
    try:
        rows = table.find_all('tr')
    except AttributeError as e:
        logging.error( 'No table rows found, exiting')
        return 1

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
