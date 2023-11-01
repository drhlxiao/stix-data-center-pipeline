import pymongo
import sqlite3
import requests
import os
from datetime import datetime
from dateutil import parser as dtparser
from stix.spice import time_utils as sdt
import sys
connect = pymongo.MongoClient('localhost')
mdb = connect["stix"]
db = mdb['eui']

def read_eui_db(filename):
    
    #last_doc = list(db.find({}).sort('start_utc',-1).limit(1))
    #start_unix =   last_doc[0]['start_unix']  if last_doc  else None
    sqlite_conn = sqlite3.connect(filename)
    sqlite_cursor = sqlite_conn.cursor()

    db.delete_many({})

    # we do a purge, it is much faster than insert many by one
    #if start_unix:
    #    start_utc=sdt.unix2utc(start_unix)
    #    sqlite_query = f'SELECT "date-avg", "imgtype", detector, wavelnth,filename,filter FROM fits_file WHERE "date-avg" > "{start_utc}";'
    #else:
    sqlite_query = f'SELECT "date-avg", imgtype, detector, wavelnth,filename,filter FROM fits_file'

    print(sqlite_query)
    sqlite_cursor.execute(sqlite_query)
    #selected_records = sqlite_cursor.fetchall()
    print('inserting data to database...')
    num=0
    new_docs=[]
    for row in sqlite_cursor:
        doc={'start_unix': sdt.utc2unix(row[0]),
                 'start_utc':row[0],
                 #'detector':row[2],
                 'imgtype':row[1],
                 'filter':row[5],
                 'wavelength':row[3],
                 'link': row[4]
                    }
        new_docs.append(doc)
        num+=1
        if num%1000==0:
            print(row[0])
            print(f'Number of records inserted: {num}')
            db.insert_many(new_docs)
            new_docs=[]
            #db.update_one({'start_utc':doc['start_utc'], 'filter':doc['filter']},{'$set': doc},upsert=True)
        

def download_file(url, local_filename='/data/EUI/metadata.db'):
    # Check if the local file exists
    has_new = True
    if os.path.exists(local_filename):
        # Get the last modified time of the local file
        local_file_time = os.path.getmtime(local_filename)
        local_file_datetime = datetime.fromtimestamp(local_file_time)

        # Send a HEAD request to get the last modified time of the remote file
        response = requests.head(url)
        remote_last_modified = response.headers.get('last-modified')
        if remote_last_modified:
            remote_file_datetime = datetime.strptime(remote_last_modified, '%a, %d %b %Y %H:%M:%S %Z')

            # Compare the last modified times
            if remote_file_datetime <= local_file_datetime:
                has_new=False
                print(f"The file '{local_filename}' is up to date. No need to download.")
                return has_new, local_filename
            else:
                print(f"The file '{local_filename}' is outdated. Downloading the updated file.")
        else:
            print("Failed to retrieve the last modified time of the remote file. Downloading anyway.")
    else:
        print("The local file does not exist. Downloading the file.")

    # Download the file
    print(f"Downloading file ...")
    response = requests.get(url)
    with open(local_filename, 'wb') as f:
        f.write(response.content)
    return has_new, local_filename
    print(f"Downloaded and saved the file as '{local_filename}'.")


def main():

    url= "https://www.sidc.be/EUI/data/metadata.db"
    filename= "/mnt/nas05/stix/EUI/metadata.db"
    has_new, fname=download_file(url, filename)
    if has_new:
        read_eui_db(fname)

if __name__ == "__main__":
    main()




