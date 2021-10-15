#export timestamps
import sys
import pymongo
import math
def to_txt(filename, server='localhost', port=27017):
    try:
        connect = pymongo.MongoClient(server, port)
        db = connect["stix"]
        collection= db['packets']
    except Exception as e:
        print('request')
        print(str(e))

    f=open(filename,'w')

    print('request data...')
    timestamps=collection.find({'_id':{'$gte':73}}, {'header.unix_time':1}).sort('header.unix_time',1)
    last_time=0
    print('writting data...')
    gaps=[]
    for i, tt  in enumerate(timestamps):
        ti=tt['header']['unix_time']
        if ti-last_time>100:
            f.write('{}\n'.format(ti))
        last_time=ti

    print('done..')
    f.close()
    print(filename)   



if __name__ == '__main__':
    port=27017
    if len(sys.argv)>=2:
        to_txt(sys.argv[1], port=port)
    else:
        print(' run filename, port')
