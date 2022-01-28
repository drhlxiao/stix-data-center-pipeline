#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pymongo
from dateutil import parser as dtparser
#fname='/home/xiaohl/FHNW/STIX/SolarFlareAnalysis/eui/fsi09.csv'
import pandas as pd
import csv
from stix.spice import datetime as sdt
import sys

connect = pymongo.MongoClient()
mdb = connect["stix"]
db = mdb['eui']


# In[22]:

def main(fname):
    with open(fname) as f:
        lines=csv.reader(f,delimiter='\t')
        #print(lines)
        rows=[row for row in lines]
        print(rows[4])
        #continue
        for row in rows[1:]:
            if not row:
                continue
            ut=sdt.utc2unix(row[0])
            try:
                link=''
                for col in row[7:]:
                    if '.fits' in col:
                        link=col
                        break
                
                doc={'start_unix': ut,
                 'start_utc':row[0],
                 'detector':row[1],
                 'filter':row[6],
                 'link':link
                    }
                print(doc)
                #docs.append(doc)
                db.update_one({'start_utc':doc['start_utc'], 'filter':doc['filter']},{'$set': doc},upsert=True)

            except Exception as e:
                print(row)
                print(e)
            #print(doc)
        


# In[14]:
if __name__=='__main__':
    print('import eui <filaneme>')
    main(sys.argv[1])



# In[ ]:




