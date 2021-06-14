#!/usr/bin/python3
import ior2dict
import time
import os
import sys
DESCRIPTION = """
A tool to validate STIX PDOR Printouts.
Usage:
    printout_validate PDOR_DIRECTORY  [PRINTOUT_DIRECTORY]
"""

import sys
import os
from pprint import pprint
import csv

ior_reader = ior2dict.IORReader()


#pprint(CSS)
def error(msg):
    print('[ERROR]:{}'.format(msg))


def info(msg):
    print('[SUCCESS]:{}'.format(msg))


def read_printouts(filename):
    occurrences = []
    with open(filename) as f:
        reader = csv.reader(f)
        occur = {}
        for row in reader:
            if row[0]:
                if row[0].isnumeric():
                    if occur:
                        occurrences.append(occur)
                    occur = {}
                    occur['name'] = row[1]
                    occur['duration'] = row[8]
                    occur['seq'] = row[15]
                    occur['parameters'] = []
            else:
                if not occur:
                    print('Parameter: {} does not have TC'.format(str(row)))
                    continue
                if row[1] == 'FIXED':
                    continue
                idb_param_name = row[1]
                if 'PIX' not in idb_param_name:  #convert parameter name to stix parameter name
                    ior_reader.get_idb_parameter_name(occur['seq'], row[1])
                param = (idb_param_name, row[3], row[5], row[7])
                #print(row)
                occur['parameters'].append(param)
        else:
            if occur:
                #last occurrence
                occurrences.append(occur)

    return occurrences


def equal(param1, param2):
    if param1[4] != param2[0]:
        return False
    value1 = param1[1]
    value2 = param2[3]
    if value1.startswith('0x'):
        value1 = int(value1, 16)
    if param2[2] == 'Hex':
        value1 = int(value1)
        value2 = int('0x' + param2[3], 16)
    if value1 == value2:
        return True
    else:
        try:
            if int(value1) == int(value2):
                return True
        except:
            error('{} not same: {} -  {}'.format(param1[4], value1, value2))
            return False


def check_parameters(params1, params2):
    dismatches = []
    for p1 in params1:
        if not p1:
            continue
        match = False
        for p2 in params2:
            if equal(p1, p2):
                match = True
                break
        if not match:
            dismatches.append(p1[4])
    return dismatches


def check_occurrence(occ, printouts):
    pro = None
    while True:
        pro = printouts.pop(0)
        if 'ZIX' in pro['name']:
            break

    if pro['name'] == occ['name']:
        paramter_mismatches = check_parameters(occ['parameters'],
                                               pro['parameters'])

        if paramter_mismatches:
            error('Mismatched parameters in {}: {}'.format(
                occ['name'], str(paramter_mismatches)))
        else:
            info('{} validated'.format(occ['name']))
    else:
        error('{} are not the same '.format(pro['name']))


def diff_pdor_printout(pdor_filename, printout_filename):
    pdors = ior_reader.parse(pdor_filename, no_seq=True)
    printouts = read_printouts(printout_filename)
    for occ in pdors['occurrences']:
        if 'ZIX' not in occ['name']:
            info('Platform command {} ignored'.format(occ['name']))
        else:
            check_occurrence(occ, printouts)


#pf='test_files/PDOR_SSTX_S001_CFTPART1_00001_F.csv'
#pdorf='test_files/PDOR_SSTX_S001_CFTPART1_00001.SOL'

#diff_pdor_printout(pdorf, pf)


def compare_requests(pdor_path, printout_path=None):
    cfiles = []
    if not printout_path:
        printout_path = pdor_path
    html_filename = os.path.join(printout_path, 'printouts_check.html')
    with open(html_filename) as html:
        for root, dirs, files in os.walk(printout_path):
            for file in files:
                if file.endswith('.csv') and file.startswith('PDOR'):
                    printout_filename = os.path.join(root, file)
                    pdor_filename = os.path.splitext(pdor_filename)[0] + '.SOL'
                    diff_pdor_printout(pdor_filename, printout_filename)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        compare_requests('./')
    elif len(sys.argv) == 2:
        compare_requests(sys.argv[1])
    else:
        compare_requests(sys.argv[1], sys.argv[2])
