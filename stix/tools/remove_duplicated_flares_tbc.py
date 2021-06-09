import pymongo
connect = pymongo.MongoClient()
db = connect["stix"]
fdb= db['flares']
cursor=fdb.find()
num=0
for f in  cursor:
     hidden=False
     time_window=600
     peak_unix=f['peak_unix_time']
     print(peak_unix)
     exists=fdb.find_one({'peak_unix_time':{'$gt':peak_unix-time_window, '$lt':peak_unix+time_window}, '_id':{'$ne':f['_id']}})
     if exists:
         hidden=True
     f['hidden']=hidden
     fdb.replace_one({'_id':f['_id']}, f)
     print(num,hidden)
     num+=1

            


    

    
    
