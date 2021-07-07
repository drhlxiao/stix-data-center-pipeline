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
from log_file import get_logfile

def insert_database(trigid, gdatetime, ra, dec,imgurl):
    #insert the information into database
    #to be implemented later
    print("Inserting into database")
    print(trigid)
    print(gdatetime)


def extract_image_url(item):
    '''Fermi web page doesn't show light curve image link directly,
    one has to extract it from another web page (a text file), whose link is provided in the table
     this function find the link, open the link and then find the line start with LC_URL: http://****
     '''
    
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
            print('An error occurred fetching %s \n %s' % (url, e.reason))   

    return  iurl 


def parse_rows(rows):
    """ Get data from rows """
    msglist= []
    nn=0

    for row in rows:
        #iterate each row in the table
        table_data = row.find_all('td')
        if len(table_data)>7:
            trigid=table_data[0].get_text()
            if nn>10:
                break
            
            if trigid.isdigit():
            #extract table cell content
                try:
                    print(trigid)
                    imgurl=extract_image_url(table_data[0])
                    gdate=table_data[1].get_text()
                    gtime=table_data[2].get_text()
                    gdatetime=datetime.strptime(gdate,'%y/%m/%d').strftime("%Y-%m-%d")+" "+gtime
                    #date time of the burst, merge date and time

                    ra=table_data[4].get_text()
                    dec=table_data[5].get_text()
                    #push to database
                    insert_database(trigid, gdatetime, ra, dec,imgurl)


                    nn=nn+1

                except:
                       continue


def main():

    url='http://gcn.gsfc.nasa.gov/fermi_grbs.html'
    try:
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

    try:
        rows = table.find_all('tr')
    except AttributeError as e:
        print('No table rows found, exiting')
        return 1

    table_data = parse_rows(rows)


if __name__ == '__main__':
    status = main()
    sys.exit(status)
