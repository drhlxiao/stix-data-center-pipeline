"""
http services for the idl imaging pipeline, 
The port is only open for local machine
Author:hualin.xiao@fhnw.ch
Date:2022-05-03
"""

import os
from datetime import datetime
from stix.core import mongo_db as mgdb
from flask import Flask,    jsonify, request
app = Flask(__name__)

mdb = mgdb.MongoDB()

@app.route('/request/imaging/task/last')
def get_last_pending_request():
    #used by idl daemon
    db=mdb.get_collection('flare_images')
    docs=list(db.find({'num_idl_calls':0} ).sort('_id',-1).limit(1))
    res={'pending':0}
    if docs:
        res=docs[0] 
        res['pending']=1
    return jsonify(res) 
@app.route('/request/imaging/task/update', methods=['POST', 'GET'])
def update_request():
    data=request.form.to_dict()
    if '_id' not in data:
        return jsonify({"error":'ID not specified'})
    _id=int(data['_id'].strip())
    fits={key: value.strip() for key,value in data.items() if key!='_id'}
    db=mdb.get_collection('flare_images')
    print("Updating: ", _id)
    db.update_one({'_id':_id},{'$set':{'fits':fits, 'num_idl_calls':1, 'idl_status':True}}, upsert=False) 
    res={'success':True}
    return jsonify(res)
if __name__=='__main__':
    app.run(port='8022')



