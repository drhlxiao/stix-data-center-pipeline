import csv
f=open('STIX_flare_list_commissioning.csv')
data = csv.DictReader(f)
import pymongo
from dateutil import parser as dtparser
connect = pymongo.MongoClient()
db = connect["stix"]
col = db['flares']

def utc2unix(utc):
    if isinstance(utc, str):
        if not utc.endswith('Z'):
            utc += 'Z'
        try:
            return dtparser.parse(utc).timestamp()
        except:
            return 0
    elif isinstance(utc, int) or isinstance(utc, float):
        return utc
    else:
        return 0
    


for i, l in enumerate(data):
    l['_id']=i
    l['start_unix_time']=utc2unix(l['date']+'T'+l['start'])
    l['duration']=-(l['start_unix_time']-utc2unix(l['date']+'T'+l['end']))
    print(l)
    col.insert_one(l)
