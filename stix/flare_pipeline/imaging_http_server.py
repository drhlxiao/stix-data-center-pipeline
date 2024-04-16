"""
http services for the idl imaging pipeline, 
The port is only open for local machine
Author:hualin.xiao@fhnw.ch
Date:2022-05-03
"""

import os
from datetime import datetime
from stix.core import mongo_db as mgdb
from flask import Flask, jsonify, request

app = Flask(__name__)

mdb = mgdb.MongoDB()


@app.route('/request/imaging/task/last')
def get_last_pending_request():
    #used by idl daemon
    db = mdb.get_collection('flare_images')
    docs = list(db.find({'num_idl_calls': 0}).sort('_id', -1).limit(1))
    res = {'pending': 0}
    if docs:
        res = docs[0]
        res['pending'] = 1
        try:
            model = str(res['idl_config']['ospex']['model'])
        except KeyError:
            model = 'auto'
        res['require_nonthermal'] = 0
        res['thermal_only'] = 0
        #auto
        if model == 'vth':
            res['thermal_only'] = 1
        elif model == 'vth+thick2':
            res['require_nonthermal'] = 1

    print(res)

    return jsonify(res)


@app.route('/request/imaging/task/update', methods=['POST', 'GET'])
def update_request():
    data = request.form.to_dict()
    if '_id' not in data:
        return jsonify({"error": 'ID not specified'})
    _id = int(data['_id'].strip())
    db = mdb.get_collection('flare_images')
    fields={}
    fits = {}
    vis={}
    images_stix_frame={}
    idl_status=data.get('idl_status','')
    error=data.get('error',None)
    idl_message='' if error is None else error
    aux_sav_fname =None
    for key, value in data.items():
        value = value.strip()
        if key in ['_id', 'idl_status', 'error'] or not value:
            continue
        if os.path.exists(value):
            if key == 'vis':
                vis['sav_filename'] = value
            elif 'image_stix' in key:
                images_stix_frame[key]= value
            elif key=='aux':
                aux_sav_fname = value
            else:
                fits[key] = value

    fields={
            'py_calls':0,
            'fits': fits,
            'image_fits_stix_frame': images_stix_frame,
            'aux_sav': aux_sav_fname,
            'vis':vis,
            'processing_date': datetime.now(),
            'num_idl_calls': 1,
            'idl_status': idl_status,
            'idl_message': idl_message
        }
    db.update_one({'_id': _id}, {
        '$set': fields
    }, upsert=False)
    res = {'success': True}
    return jsonify(res)


def run_server(port=8022):
    app.run(port=port)


if __name__ == '__main__':
    run_server()
