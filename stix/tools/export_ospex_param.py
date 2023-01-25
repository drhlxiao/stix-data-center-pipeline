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
mdb = db.MongoDB()
flare_image_db = mdb.get_collection('flare_images')

from ROOT import TFile, TTree

print('function definition')

def export_meta():
    fname='ospex_meta.root'
    print('creating root file:', fname)
            
    f=TFile(fname,"recreate")
    t=TTree("tree","tree")
    is_set=False
    branch_names={}
    nmax=10
    emodel=np.array([0], dtype=np.int32)

    t.Branch('emodel',emodel,'emodel/I')

    params=np.array([0]*9, dtype=np.float64)

    t.Branch('params',params,'params[9]/D')
    param_err=np.array([0]*9, dtype=np.float64)

    t.Branch('param_err',param_err,'param_err[9]/D')

    chi2=np.array([0], dtype=np.float64)
    t.Branch('chi2',chi2,'chi2/D')
    print('looping over docs')

    for doc in flare_image_db.find({'ospex_meta':{'$exists':True}}).sort('_id', -1):
        print(doc['_id'])
        flare_aux=doc['aux']
        ospex=doc['ospex_meta']
        emodel[0]=0
        if 'auto' in ospex['fitfun']:
            emodel[0]=2
        elif 'thick2' in ospex['fitfun']:
            emodel[0]=1 

        for i,param in enumerate(ospex['params']):
            try:
                if float(param)==param:
                    #is number
                    params[i]=param
            except:
                pass
        for key, value in flare_aux.items():
            try:
                value=float(value)
            except (TypeError, ValueError):
                continue

            if not is_set and key not in branch_names:
                branch_names[key]=np.array([0],np.float64)
                t.Branch(key,branch_names[key],f'{key}/D')
            branch_names[key][0]=value

        chi2[0]=ospex['chi2']
        is_set = True
        
        t.Fill()
    t.Write()
    f.Close()

if __name__=='__main__':
    export_meta()
