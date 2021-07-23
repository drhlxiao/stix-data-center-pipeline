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
from datetime import datetime,timedelta
from  grb_db import *
from log_file import get_logfile

import logging





def insert_database(trigid, gdatetime, ra, dec,url):
    
    msg=""
    if not grb_exists('GBM', trigid, gdatetime):
        logging.info("inserting...")
        ins="insert into burst_events(name, start_utc, stop_utc, burst_utc, right_ascension, declination, url, comments, author) values('GBM %s', '%s','%s','%s', %f,%f,'%s','Also detected by Fermi GBM (trigger#%s)','robot')"%(trigid, gdatetime,gdatetime,gdatetime,float(ra),float(dec),url,trigid)
        excute_sql(ins)
        msg="GBM %s    %s   %f   %f"%(trigid, gdatetime, float(ra),float(dec))
    else:
        logging.info("Already existed fermi trigger:")
        logging.info("                   :%s \t %s \t %s \t %s"%(trigid, gdatetime,ra, dec))

    
    return msg


def parse_rows(rows):
    """ Get data from rows """
    nn=0

    for row in rows:
        table_data = row.find_all('td')
        if True:
            try:
                trigid=table_data[0].get_text()
                logging.info(trigid)

                imgurl=''
                gdate=table_data[1].get_text()
                gtime=table_data[2].get_text()
                da=datetime.strptime(gdate,'%y/%m/%d')
                #if da < datetime.now()-timedelta(days=360):
                #    continue
                gdatetime=da.strftime("%Y-%m-%d")+" "+gtime
                print(gdatetime)

                ra=table_data[4].get_text()
                dec=table_data[5].get_text()
                msg=insert_database(trigid, gdatetime, ra, dec,imgurl)

                #ut=dtparser.parse(gdatetime).timestamp()
                #start=int(ut)-900
                #url=f'https://pub023.cs.technik.fhnw.ch/view/plot/lightcurves?start={start}&span=1800'
                #urls.append([url,row[1],row[0]])
                nn=nn+1

            except Exception as e:
                print(e)


def main():
    # Get arguments
    # Make soup
    url='http://gcn.gsfc.nasa.gov/fermi_grbs.html'
    try:
        logging.info("opening %s"%url)
        resp = urlopen(url)
    except URLError as e:
        logging.info('An error occured fetching %s \n %s' % (url, e.reason) )
        return 1
    soup = BeautifulSoup(resp.read(),"html.parser")
    # Get table
    try:
        table = soup.find('table')
    except AttributeError as e:
        logging.info('No tables found, exiting')
        return 1

    # Get rows
    try:
        rows = table.find_all('tr')
    except AttributeError as e:
        logging.info('No table rows found, exiting')
        return 1
    print("RUNS:",len(rows))

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
