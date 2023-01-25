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
from stix.core import mongo_db as db
from stix.core import logger
import numpy as np

mdb = db.MongoDB()
flare_image_db = mdb.get_collection('flare_images')
logger = logger.get_logger()





def export_meta():
    fname='image_ospex_meta.root'
    print('creating root file:', fname)
            
    from ROOT import TFile, TTree
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

    
    nmax=10
    npeaks=np.array([0], np.int32)
    t.Branch('npeaks',npeaks,'npeaks/I')

    is_image_ok=np.array([0], np.int32)
    t.Branch('is_image_ok',is_image_ok,'is_image_ok/I')


    area=np.array([0]*nmax, np.float64)

    t.Branch('area',area,'area[10]/D')
    peak_x=np.array([0]*nmax, np.float64)

    t.Branch('peak_x',peak_x,'peak_x[10]/D')

    peak_y=np.array([0]*nmax,np.float64)
    t.Branch('peak_y',peak_y,'peak_y[10]/D')

    rsun=np.array([0],np.float64)
    t.Branch('rsun',rsun,'rsun[10]/D')


    print('looping over docs')

    for doc in flare_image_db.find({'ospex_meta':{'$exists':True},  'image_meta.flare_aux':{'$exists':True} }).sort('_id', -1):
        flare_aux=doc['image_meta']['flare_aux']
        print(doc['_id'])
        is_image_ok[0]=0
        sources=doc['image_meta']['sources']
        ospex=doc['ospex_meta']
        emodel[0]=0
        rsun[0]=doc['aux']['rsun']
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
        npeaks[0]=len(sources)
        for k, source in enumerate(sources):
            if k>=nmax:
                break
            try:
                area[k]=source['contour']['area']
                peak_x[k]=source['peaks'][0][0]
                peak_y[k]=source['peaks'][0][1]
                is_image_ok[0]=1
            except: 
                is_image_ok[0]=0
                pass


        chi2[0]=ospex['chi2']
        is_set = True
        
        t.Fill()
    t.Write()
    f.Close()

if __name__=='__main__':
    export_meta()
