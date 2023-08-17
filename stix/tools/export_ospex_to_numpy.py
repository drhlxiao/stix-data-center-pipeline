#!/usr/bin/env python3  
# -*- coding: utf-8 -*- 
#----------------------------------------------------------------------------
# Created By  : Hualin Xiao (hualin.xiao)
# Created Date: Sept 5, 2022
# version =1.0
"""
    Various tools to manipulate stix images
    python version must be greater than  3.9
"""
import os
import sys
import numpy as np
sys.path.append('.')

from stix.core import mongo_db as db
mdb = db.MongoDB(port=9000)
flare_image_db = mdb.get_collection('flare_images')

import json

print('function definition')

def export_meta():
    fname='ospex_meta.root'
    print('creating root file:', fname)
            
    is_set=False
    branch_names={}
    nmax=10

    docs=flare_image_db.find({'ospex_meta':{'$exists':True}}, {'ospex_meta':1, 'start_unix':1, 'end_unix':1,
        'energy_range':1, 'total_counts':1, '_id':1, 'fits.ospex_results':1, 'unique_id':1
        })
        
    with open('flare_db.json', 'w') as f:
        data=list(docs)
        json.dump(data, f)

if __name__=='__main__':
    export_meta()
