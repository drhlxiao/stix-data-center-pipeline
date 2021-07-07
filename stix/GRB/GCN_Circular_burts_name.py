#!/usr/bin/python
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



'''
def grb_exists(burst_time):
    sql='select * from burst_events where abs(timestampdiff(second, burst_utc,"'+burst_time+'"))<10 '
    res=get_rowcount_sql(sql)
    if res>0:
        return True
    else:
        return False
'''


def insert_database(trigid, gdatetime, ra, dec,url):
    
    msg=""
    if not grb_exists('Fermi GBM (trigger#', trigid, gdatetime):
        print("inserting...")
        ins="insert into burst_events(name, start_utc, stop_utc, burst_utc, right_ascension, declination, url, comments, author) values('GBM %s', DATE_SUB('%s',INTERVAL '30' SECOND), DATE_ADD('%s',INTERVAL '180' SECOND), '%s', %f,%f,'%s','Also detected by Fermi GBM (trigger#%s)','robot')"%(trigid, gdatetime,gdatetime,gdatetime,float(ra),float(dec),url,trigid)
        excute_sql(ins)
        msg="GBM %s    %s   %f   %f"%(trigid, gdatetime, float(ra),float(dec))
    else:
        print("Already existed fermi trigger:")
        print("                   :%s \t %s \t %s \t %s"%(trigid, gdatetime,ra, dec))

    
    return msg

def extract_image_url(item):
    iurl=""
    for link in item.find_all('a'):
        url="http://gcn.gsfc.nasa.gov/"+link.get('href')
        try:
            resp = urlopen(url)
            for line in resp:
                urls = re.findall('LC_URL:\s*http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*\(\),]|(?:%[0-9a-fA-F][0-9a-fA-F]))+', line)
                if urls:
                    iurl=re.sub('LC_URL:','',urls[0]).strip()
                    return iurl

        except URLError as e:
            print('An error occured fetching %s \n %s' % (url, e.reason))   

    return  iurl 


def parse_rows(rows):
    """ Get data from rows """
    nn=0

    for row in rows:
        table_data = row.find_all('td')
        if len(table_data)>7:
            trigid=table_data[0].get_text()
            if nn>10:
                break
            
            if trigid.isdigit():
                try:
                    print(trigid)
                    imgurl=extract_image_url(table_data[0])
                    gdate=table_data[1].get_text()
                    gtime=table_data[2].get_text()
                    gdatetime=datetime.strptime(gdate,'%y/%m/%d').strftime("%Y-%m-%d")+" "+gtime
                    ra=table_data[4].get_text()
                    dec=table_data[5].get_text()
                    msg=insert_database(trigid, gdatetime, ra, dec,imgurl)

                    nn=nn+1

                except:
                       continue


def main():
    # Get arguments
    # Make soup
    url='http://gcn.gsfc.nasa.gov/fermi_grbs.html'
    try:
        print("opening %s"%url)
        resp = urlopen(url)
    except URLError as e:
        print('An error occured fetching %s \n %s' % (url, e.reason))   
        return 1
    soup = BeautifulSoup(resp.read(),"html.parser")
    # Get table
    try:
        table = soup.find('table')
    except AttributeError as e:
        print('No tables found, exiting')
        return 1

    # Get rows
    try:
        rows = table.find_all('tr')
    except AttributeError as e:
        print('No table rows found, exiting')
        return 1

    table_data = parse_rows(rows)


if __name__ == '__main__':
    status = main()
    sys.exit(status)
