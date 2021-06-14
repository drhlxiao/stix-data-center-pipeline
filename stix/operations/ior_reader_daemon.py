#!/usr/bin/python3
# ior reader daemon
import ior2mongo
import time
import os
import sys

PDOR_DIRECTORIES = [
    '/home/xiaohl/FHNW/STIX/NECP/2020-04-22-NECP_IX-2_Rev2/PDORs/'
]
MAX_SIZE = 10 * 1024 * 1024
#1MB maxima

file_md5 = {}
MDB = ior2mongo.IORMongoDB()


def find_files(path):
    cfiles = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.SOL') and (file.startswith('PDOR')
                                          or file.startswith('IOR')):
                filename = os.path.join(root, file)
                if os.path.getsize(filename) <= MAX_SIZE:
                    md5 = ior2mongo.compute_md5(filename)
                    if md5 in file_md5.values():
                        #not to process again
                        continue
                    cfiles.append(filename)
                    file_md5[filename] = md5

                else:
                    print("Big file ignored: {}".format(filename))
    return cfiles


def update_ior_description(paths=PDOR_DIRECTORIES):
    """ read description from man.csv
        not used any more
    """

    for folder in paths:
        csv_filename = os.path.join(folder, 'man.csv')
        if not os.path.exists(csv_filename):
            continue
        md5 = ior2mongo.compute_md5(csv_filename)
        print('find man.csv in folder ', folder)
        if md5 in file_md5.values():
            print('man.csv processed.')
            continue
        print('processing manifest file...')
        with open(csv_filename) as f:
            lines = f.readlines()
            for line in lines:
                columns = line.split(',')
                if len(columns) == 3:
                    print('updating :')
                    print(columns)
                    MDB.update_ior_info(columns[0],
                                        description=columns[2],
                                        phase=columns[1])
        file_md5[csv_filename] = md5


def scan_ior(paths=PDOR_DIRECTORIES):
    for path in paths:
        print("Checking directory:", path)
        filelist = find_files(path)
        for filename in filelist:
            MDB.insert(filename)
        else:
            print('No new PDOR found')


def main_loop():
    while True:
        scan_ior()
        print('wait for 60 seconds ...')
        time.sleep(60)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main_loop()
    else:
        path = sys.argv[1]
        answer = input("Import PDORs in {} ?".format(path))
        if answer in ['Y', 'y']:
            scan_ior([path])
            #update_ior_description([path])
