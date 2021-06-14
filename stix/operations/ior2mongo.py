#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : PDOR2MongoDB.py
# @description  : import IORs and dump them into MongoDB
# @author       : Hualin Xiao
# @date         : Jan. 22, 2020
#import json
import pprint
import sys
import os
import pymongo
import hashlib
from datetime import datetime
from dateutil import parser as dtparser
import ior2dict
MIB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        '../data/MIB')


def get_now(dtype='unix'):
    utc_iso = datetime.utcnow().isoformat() + 'Z'
    return dtparser.parse(utc_iso).timestamp()


def compute_md5(filename):
    file_hash = hashlib.md5()
    with open(filename, "rb") as f:
        while True:
            chunk = f.read(8192)
            if not chunk:
                break
            file_hash.update(chunk)
    return file_hash.hexdigest()


class IORMongoDB(object):
    def __init__(self,
                 server='localhost',
                 port=27017,
                 user='',
                 pwd='',
                 MIB=MIB_PATH):
        self.filename = None
        self.packets = []
        self.db = None
        self.collection_IOR = None
        self.current_id = 0
        self.ior_reader = ior2dict.IORReader(MIB)
        try:
            if server == 'localhost' and user == '' and pwd == '':
                self.connect = pymongo.MongoClient(server, port)
            else:
                self.connect = pymongo.MongoClient(server,
                                                   port,
                                                   username=user,
                                                   password=pwd,
                                                   authSource='stix')
            self.db = self.connect["stix"]
            self.collection_IOR = self.db['iors']
        except Exception as e:
            print('can not connect to mongodb')

    def update_ior_info(self, ior_filename, description, phase=''):
        filename = os.path.basename(ior_filename)
        if self.collection_IOR:
            ior = self.collection_IOR.find_one({'filename': filename})
            if ior:
                ior['description'] = description
                ior['phase'] = phase
                self.collection_IOR.save(ior)

    def is_connected(self):
        if self.db:
            return True
        else:
            return False

    def is_processed(self, md5):
        try:
            #print(md5)
            run = self.collection_IOR.find_one({'md5': md5})
            if run:
                return True
        except Exception as e:
            print(str(e))
        return False

    def insert(self, in_filename):
        if not self.db:
            print('MongoDB not connected')
            return
        with open(in_filename) as fin:
            md5 = compute_md5(in_filename)
            filename = os.path.basename(in_filename)
            abspath = os.path.abspath(in_filename)
            path = os.path.dirname(abspath)
            if self.is_processed(md5):
                print('The IOR {} is already in the database'.format(
                    in_filename))
                return
            try:
                self.current_id = self.collection_IOR.find().sort(
                    '_id', -1).limit(1)[0]['_id'] + 1
            except IndexError:
                self.current_id = 0
            request = {
                'filename': filename,
                'path': path,
                'md5': md5,
                'description': '',
                'hidden': False,
                'phase': '',
                'processing_time': get_now(),
                'status': 1,
                'log': ''
            }

            data = {}
            try:
                data = self.ior_reader.parse(in_filename)
                request.update(data)
            except Exception as e:
                print(str(e))
                request['log'] = str(e)

            request['_id'] = self.current_id
            self.collection_IOR.insert_one(request)
            print('Insert successfully')
            print('Entry ID:{}'.format(self.current_id))


if __name__ == '__main__':
    if len(sys.argv) == 2:
        mdb = IORMongoDB()
        mdb.insert(sys.argv[1])
    else:
        print('IOR2MongoDB <PDOR.SOL>')
