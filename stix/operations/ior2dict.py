#!/usr/bin/python3
# A script to extract information from PDOR files
#@Author: Hualin Xiao
#@Date:   Feb. 11, 2020

import sys
import os
from pprint import pprint
import re
import xml.etree.ElementTree as ET
from datetime import datetime
from dateutil import parser as dtparser
from copy import deepcopy
MIB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        '../data/MIB')


def to_epoch(time_str, dformat='%Y-%jT%H:%M:%SZ'):
    #time_str format:  2049-077T18:59:59Z
    # create 1,1,1970 in same timezone as d1
    ts = 0
    try:
        d1 = datetime.strptime(time_str, dformat)
        d2 = datetime(1970, 1, 1, tzinfo=d1.tzinfo)
        time_delta = d1 - d2
        ts = int(time_delta.total_seconds())
    except Exception as e:
        print(str(e))
        pass

    return ts


def remove_namespace(doc, namespace):
    """Remove namespace in the passed document in place."""
    ns = u'{%s}' % namespace
    nsl = len(ns)
    for elem in doc.getiterator():
        if elem.tag.startswith(ns):
            elem.tag = elem.tag[nsl:]


def unix2utc(ts):
    return datetime.utcfromtimestamp(ts).isoformat()


def get_aspect_uid(utc):
    if not utc.endswith('Z'):
        utc += 'Z'
    level = 5
    version = 0
    start_datetime = dtparser.parse(utc)
    year = start_datetime.year
    month = start_datetime.month
    day = start_datetime.day
    hour = start_datetime.hour
    minute = start_datetime.minute
    request_id = (year & 0xF) << 28
    request_id |= (month & 0xF) << 24
    request_id |= (day & 0x1F) << 19
    request_id |= (hour & 0x1F) << 14
    request_id |= (minute & 0x3F) << 8
    request_id |= (level & 0xF) << 4
    request_id += version
    return request_id


class IORReader(object):
    def __init__(self, mib_path=MIB_PATH):

        self.table_names = ['ccf', 'csf', 'sdf', 'cdf']
        self.mib = {}
        self.read_MIB(mib_path)
        self.read_css(mib_path)
        self.tc_time_offset = 0

    def read_MIB(self, mib_path):
        for table in self.table_names:
            filename = os.path.join(mib_path, table + '.dat')
            self.mib[table] = []
            with open(filename) as f:
                lines = f.readlines()
                for line in lines:
                    self.mib[table].append(line.split('\t'))

    def get_idb_parameter_name(self, sequence, name):
        if 'PIX' in name:
            return name
        for line in self.mib['sdf']:
            if sequence == line[0] and name == line[7]:
                return line[4]

        return ''

    def get_telecommand_parameters(self, telecommand):
        return [x[6] for x in self.mib['cdf'] if telecommand == x[0] and x[6]]

    def get_description(self, command):
        if 'ZIX' in command:
            for tc in self.mib['ccf']:
                if tc[0] == command:
                    return tc[1:3]
        if 'AIX' in command:
            for tc in self.mib['csf']:
                if tc[0] == command:
                    return tc[1:3]
        return ['', '']

    def find_parameter(self, name, occ_parameters):
        if not occ_parameters:
            return []
        for occ_param in occ_parameters:
            if name == occ_param[0] or name == occ_param[-1]:
                return occ_param
        return []

    def get_action_time(self, start, tc):
        #print(start, tc)
        time_type = tc[3]  # A: absolute, R: relative
        action_time = tc[4]
        action_time_seconds = tc[5]
        is_absolute = re.findall(r'\d{4}-\d{2}-\d{2}', start)
        if is_absolute:
            unix_time = dtparser.parse(start).timestamp()
            if time_type == 'A':
                return start
            return datetime.utcfromtimestamp(unix_time +
                                             action_time_seconds).isoformat(
                                                 timespec='milliseconds')
        else:
            if time_type == 'R':
                return action_time
            else:
                return start

        return ''

    def seq2tc(self, occurrences):
        #convert all sequences to occurrences
        new_occurrences = []
        for occ in occurrences:
            if 'ZIX' in occ['name']:
                new_occurrences.append(deepcopy(occ))
            else:
                telecommands = self.get_telecommands(occ['name'])
                start = occ['actionTime']
                self.tc_time_offset = 0
                for tc in telecommands:
                    name = tc[0]
                    parameters = []
                    param_list = self.get_telecommand_parameters(name)
                    if param_list:
                        parameters = [
                            self.find_parameter(name, occ['parameters'])
                            for name in param_list
                        ]
                    new_occurrences.append({
                        'name':
                        name,
                        'desc':
                        tc[1],
                        'type':
                        occ['type'],
                        'actionTime':
                        self.get_action_time(start, tc),
                        'uniqueID':
                        occ['uniqueID'],
                        'sequence':
                        occ['name'],
                        'parameters':
                        parameters
                    })
        return new_occurrences

    def read_css(self, mib_path):
        output = {}
        filename = os.path.join(mib_path, 'css.dat')
        with open(filename) as f:
            lines = f.readlines()
            last_seq = ''
            for l in lines:
                line = l.split('\t')
                if line[0] not in output:
                    output[line[0]] = []

                timestamp = line[8].strip().split('.')
                delta_time = ''
                if len(timestamp) == 3:
                    delta_time = 3600 * int(timestamp[0]) + 60 * int(
                        timestamp[1]) + int(timestamp[2])
                output[line[0]].append(
                    (line[4], line[1], line[2], line[7], line[8], delta_time))
                #
        self.css = output

    def get_telecommands(self, name):
        #TC name, description, value, 'absolute or relative', deleta time, delta time in seconds
        try:
            return self.css[name]
        except KeyError:
            return []

    def parse(self, filename, no_seq=False):
        #no_seq: No sequence, convert all sequences to TCs
        request = {}
        tree = ET.parse(filename)

        root = tree.getroot()
        remove_namespace(root, 'soc.solarorbiter.org')
        #remove namespace for IOR
        req = root[0]

        header = req.find('header')
        occurrence_list = req.find('occurrenceList')
        if header.attrib['type'] not in ['PDOR', 'IOR']:
            print('Not support type')
            return None
        request['type'] = header.attrib['type']
        gen_time = header.find('genTime').text
        request['genTime'] = gen_time
        validity = header.find('validityRange')
        start_time = validity.find('startTime').text
        stop_time = validity.find('stopTime').text
        request['stopTime'] = stop_time
        request['startTime'] = start_time
        request['startUnix'] = to_epoch(start_time)
        request['stopUnix'] = to_epoch(stop_time)
        request['genUnix'] = to_epoch(gen_time)
        request['count'] = occurrence_list.attrib['count']
        request['author'] = occurrence_list.attrib['author']
        request['data_request_unique_ids'] = []
        occurrences = []
        for x in occurrence_list:
            name = x.attrib['name']
            desc = self.get_description(name)
            execution_time = None
            action_time = 'asap'
            if x.find('releaseTime'):
                try:
                    action_time = x.find('releaseTime').find('actionTime').text
                except:
                    action_time = 'Error:-1'
            elif x.find('executionTime'):
                try:
                    execution_time = x.find('executionTime').find(
                        'actionTime').text
                except:
                    action_time = 'Error:-2'
                if execution_time:
                    dformat = '%Y-%jT%H:%M:%S.%fZ'
                    if request['type'] == 'IOR':
                        dformat = '%Y-%jT%H:%M:%SZ'
                    try:
                        unix_time = to_epoch(execution_time, dformat)
                        action_time = unix2utc(unix_time)
                    except:
                        action_time = 'Error:-3'

            occ = {
                'name': name,
                'type': x.tag,
                'desc': desc,
                'uniqueID': x.find('uniqueID').text,
                'actionTime': action_time,
                'parameters': []
            }

            parameter_list = x.find('parameterList')
            if parameter_list:
                for p in parameter_list:
                    value = p.find('value')
                    param_value = value.text
                    if 'radix' in value.attrib:
                        if value.attrib['radix'] == 'Hexadecimal':
                            param_value = '0x' + param_value
                    parameter = [
                        p.attrib['name'], param_value,
                        p.find('description').text,
                        value.attrib['representation'],
                        self.get_idb_parameter_name(name, p.attrib['name'])
                    ]
                    if parameter[0] in ['PIX00076', 'XF417B01', 'XF417A01']:
                        request['data_request_unique_ids'].append(
                            int(parameter[1]))

                    occ['parameters'].append(parameter)

            occurrences.append(occ)
        if no_seq:
            #convert all sequences to TCs
            request['occurrences'] = self.seq2tc(occurrences)
        else:
            request['occurrences'] = occurrences

        return request


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('A script to parse PDOR')
        print('Usage: python parse_PDOR <PDOR filename>')
    else:
        reader = IORReader()
        ret = reader.parse(sys.argv[1])
        pprint(ret)
