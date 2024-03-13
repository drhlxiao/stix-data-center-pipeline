#find bulk science data that don't have fits file created    
# hualin.xiao created on Feb. 11, 2022
import sys
import os
import json
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib import pyplot as plt
from stix.core import mongo_db as db
from stix.spice import time_utils as sdt
from stix.fits import fits_creator
import matplotlib.ticker as mticker
from pathlib import Path
import pymongo
mdb = db.MongoDB()
connect = pymongo.MongoClient()
db = connect["stix"]
col_bsd=db['bsd']
col_fits=db['fits']
pipeline = [
    {
        "$lookup": {
            "from": "fits",
            "let": {"bsdUniqueId": "$unique_id"},
            "pipeline": [
                {
                    "$match": {
                        "$expr": {
                            "$and": [
                                {"$eq": ["$request_id", "$$bsdUniqueId"]},
                                {"$eq": ["$level", "xray-cpd"]}
                            ]
                        }
                    }
                }
            ],
            "as": "fits"
        }
    },
    {
        "$match": {
            "fits": {"$size": 0}
        }
    },
    {
        "$project": {
            "_id": 1,  # Exclude _id field
            "unique_id": 1,  # Include unique_id field
            "start_unix_time":1
        }
    },
    {"$sort": {"_id": -1}}
]

# Execute the aggregation pipeline
docs = col_bsd.aggregate(pipeline)
f=open('L1_missing.csv','w')
f.write("Unique_ID,BSD_DB_Entry_ID, Data Start UTC\n")
num=0
for doc in docs:
    line=f'{doc["unique_id"]},{doc["_id"]},{sdt.unix2utc(doc["start_unix_time"])}\n'
    print(line)
    num+=1
    f.write(line)

print("Total files:",num)

