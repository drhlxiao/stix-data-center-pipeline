#!/usr/bin/python
#Fetch external triggers
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





def insert_database(utc,url):
    logging.info( 'inserting events: '+utc)
    dd=re.split(':|-| ', utc)
    trigid="".join(dd)
    
    if not grb_exists('SPIACS', trigid,utc):
        logging.info( "inserting...")
        logging.info( 'trigid:'+trigid)
        gdatetime=utc
    
        ins="insert into burst_events(name, start_utc, stop_utc, burst_utc, right_ascension, declination, url, comments, author) values('SPIACS %s', '%s', '%s', '%s', %f,%f,'%s','Also detected by SPI ACS','robot')"%(trigid, gdatetime,gdatetime,gdatetime,-1,-1,url)
        excute_sql(ins)
        return False 
    else:
        logging.info( "Already existed ")
        return True

    



def main():
    # Get arguments
    # Make soup

    ds=time.gmtime()
    #url='http://gcn.gsfc.nasa.gov/fermi_grbs.html'
    mon="%s"%ds.tm_mon
    mon2=mon.rjust(2,'0')
    months="month=%s-%s"%(ds.tm_year,mon2)
    
    url='http://www.isdc.unige.ch/integral/ibas/cgi-bin/ibas_acs_web.cgi?%s&showall=on'%months


    try:
        logging.info( "opening %s"%url)
        resp = urlopen(url)

    except URLError as e:
        logging.error( 'An error occurred while fetching %s \n %s' % (url, e.reason))
        return 1
    soup = BeautifulSoup(resp.read(),"html.parser")
    # Get table
    try:
        table= soup.find('table',{"class":"isdctable"})
    except AttributeError as e:
        logging.info( 'No tables found, exiting')
        return 1

    # Get rows
    try:
        rows = table.find_all('tr')

    except AttributeError as e:
        logging.info( 'No table rows found, exiting')
        return 1
    # Get rows
    try:
        logging.info( 'parsing rows')
        parse_row(rows)
    except AttributeError as e:
        logging.info( 'error when parsing rows')
        return 1

def format_comment(stype,comment):
    content=""
    if stype.rstrip() is not "": 
        content=content+stype.rstrip()
    if content!="":
        content=content+','
    if comment.rstrip() is not "":
        content=content+comment.rstrip()

    return content
        

def parse_row(rows):
    nn=0

    for row in rows:
        table_data = row.find_all('td')
        logging.info( '=======================')
        if len(table_data)==9:
            if nn>10:
                break
            nn=nn+1

            utc=table_data[0].get_text()
            spi_type=table_data[2].get_text()
            comment=table_data[7].get_text()
            spi_edit=table_data[8]
            #lcurl=extract_image_url(spi_edit)
            lcurl=extract_image_unzoom_url(spi_edit)
            butc=re.findall(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}",utc)
            spi_comment=format_comment(spi_type,comment)

            if len(butc)>0:
                burst_time=str(butc[0])
                logging.info( 'checking burst: \n %s, %s, %s'%(burst_time, spi_type.rstrip(), lcurl.rstrip()))

                logging.info( 'updating /inserting database ...')
                if 'Spike' in spi_type:
                    logging.info( 'this is a spike, updating database')
                    set_not_detected(burst_time)

                else:
                    exist=insert_database(burst_time,lcurl)
                    if exist:
                        remark(burst_time,spi_comment)

                    
def set_not_detected(burst_time):
    sql='update burst_events set detected=0, comments="SPI Spike, set to not be detected by robot" where name like "%%SPI%%" and burst_utc = "%s"'%burst_time
    excute_sql(sql)
def remark(burst_time, cmt):
    sql='update burst_events set comments="%s, Also detected by SPI ACS at %s, remarked by robot" where name like "%%SPI%%" and burst_utc="%s"'%(cmt,burst_time,burst_time)
    excute_sql(sql)


def extract_image_url(item):
    iurl=""
    for link in item.find_all('a'):
        prefix="http://www.isdc.unige.ch/"
        url=prefix+link.get('href')
        try:
            resp = urlopen(url)
            soup = BeautifulSoup(resp.read(),"html.parser")
            imglist=[x['src'] for x in soup.findAll('img')]
            for lc in imglist:
                if 'triggers' in lc:
                    return prefix+lc

        except URLError as e:
            logging.error( 'An error occured fetching %s \n %s' % (url, e.reason)   )

    return  iurl 

def extract_image_unzoom_url(item):
    iurl=""
    for link in item.find_all('a'):
        prefix="http://www.isdc.unige.ch/"
        url=prefix+link.get('href')
        try:
            resp = urlopen(url)
            soup = BeautifulSoup(resp.read(),"html.parser")
            imglist=[x['href'] for x in soup.findAll('a',href=True)]
            for lc in imglist:
                if 'triggers' in lc and 'png' in lc:
                    logging.info( 'LC url:'+prefix+lc)
                    return prefix+lc

        except URLError as e:
            logging.error( 'An error occurred while fetching %s \n %s' % (url, e.reason))





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



    main()
