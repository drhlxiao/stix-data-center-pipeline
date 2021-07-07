#!/usr/bin/python

#Fetch fermi triggers
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
from log_file import get_logfile


import logging



def insert_database(trigid, gdatetime, ra, dec,url):
    msg=""
    if not grb_exists('CALETGBM', trigid, gdatetime):
        logging.info("inserting...")
        logging.info("inserting...")
        ins="insert into burst_events(name, start_utc, stop_utc, burst_utc, right_ascension, declination, url, comments, author) values('CALETGBM %s', '%s','%s', '%s', %f,%f,'%s','Also detected by CALET_GBM (trigger#%s)','robot')"%(trigid, gdatetime,gdatetime,gdatetime,float(ra),float(dec),url,trigid)
        excute_sql(ins)
        msg="GBM %s    %s   %f   %f"%(trigid, gdatetime, float(ra),float(dec))
    else:
        logging.info("Already existed trigger:")
        logging.info("                   :%s \t %s \t %s \t %s"%(trigid, gdatetime,ra, dec))

    
    return msg

def extract_image_url(item):
    iurl=""
    for link in item.find_all('a'):
        url="http://gcn.gsfc.nasa.gov/"+link.get('href')
        iurl=url
    return  iurl 


def parse_rows(rows):
    """ Get data from rows """
    nn=0

    for row in rows:
        table_data = row.find_all('td')
        if len(table_data)>7:
            trigid=table_data[0].get_text()
            if nn>20:
                break
            nn=nn+1
            
            if trigid.isdigit():
                try:
                    logging.info(trigid)
                    imgurl=extract_image_url(table_data[0])
                    gdate=table_data[1].get_text()
                    gtime=table_data[2].get_text()
                    gdatetime=datetime.strptime(gdate,'%y/%m/%d').strftime("%Y-%m-%d")+" "+gtime
                    ra=-1
                    dec=-1
                    msg=insert_database(trigid, gdatetime, ra, dec,imgurl)


                except:
                       continue


def main():
    # Get arguments
    # Make soup
    url='https://gcn.gsfc.nasa.gov/calet_triggers.html'
    try:
        logging.info('opening:%s'%url)
        resp = urlopen(url)
    except URLError as e:
        logging.error('An error occurred when fetching:%s'%url)
        return 1
    soup = BeautifulSoup(resp.read(),"html.parser")
    # Get table
    try:
        table = soup.find('table')
    except AttributeError as e:
        logging.error('No tables found')
        return 1

    # Get rows
    try:
        rows = table.find_all('tr')
    except AttributeError as e:
        logging.error('No tables found,exiting')
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
