import pymongo

connect = pymongo.MongoClient()
db = connect['stix']['auto_flares']
