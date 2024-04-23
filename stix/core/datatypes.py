#!/usr/bin/env python3

# -*- encoding: utf-8 -*-
# @title        : datatypes.py
# @date         : Oct. 28, 2019
# @description:
#               definitions of structures of decoded parameters and packets
import os
import pprint
import math
import copy
from stix.core import logger
from stix.core import idb
STIX_IDB = idb.idb()
logger = logger.get_logger()

QL_SPID_MAP = {
        54118: 'lc',
        54119:'bkg',
        54120:'qlspec',
        54121:'var',
        54122:'flare',
        }

SPID_MAP = {
        # Bulk Science
        54114: 'xray-rpd',
        54115: 'xray-cpd',
        54116: 'xray-scpd',
        54117: 'xray-vis',
        54143: 'xray-spec',
        54125: 'aspect',
        # Quick look
        54118: 'ql-lc',
        54119: 'ql-bkg',
        54120: 'ql-spec',
        54121: 'ql-var',
        54122: 'ql-ffl',
        54123: 'ql-tm',
        54124: 'ql-cal',
        # House keeping
        54101: 'hk-mini',
        54102: 'hk-maxi'
        }
LOW_LATENCY_TYPES={
        54102: 'hk-maxi',
        54118: 'ql-lc',
        54119: 'ql-bkg',
        54120: 'ql-spec',
        54121: 'ql-var',
        54122: 'ql-ffl'}

SCI_REPORT_SPIDS=[
        54114,
        54115,
        54116,
        54117,
        54143,
        54125
        ]


def copy_object(x, deep_copy):
    if deep_copy:
        return copy.deepcopy(x)
    return x


def any_list_contains(parents, child):
    for parent in parents:
        if list_contains(parent, child):
            return True
    return False


def list_contains(parent, child):
    if not child or not parent:
        return False

    if len(child) > len(parent):
        return False
    for i, e in enumerate(child):
        if parent[i] != child[i]:
            return False
    return True


class Parameter(object):

    _fields = ('name', 'raw', 'eng', 'children')

    def __init__(self, parameter=None):
        self._name = ''
        self._raw = ''
        self._eng = ''
        self._aux = ''
        self._children = []
        if isinstance(parameter, (tuple, list)):
            if len(parameter) ==4:
                self._name, self._raw, self._eng, self._children = parameter
            elif len(parameter) == 5:
                self._name, self._raw, self._eng, self._children, self._aux = parameter
                #add one parameter to keep errors in 202404
            else:
                logger.warn(
                        "Invalid parameter type, parameter value:{}".format(
                            str(parameter)))

    def copy(self, name, raw, eng, children, error=None):
        self._name = name
        self._raw = raw
        self._eng = eng
        self._children = children
        self._aux=error

    def get_raw_int(self):
        try:
            return int(self._raw)
        except (TypeError, IndexError, ValueError):
            return None

    def as_tuple(self):
        return (self._name, self._raw, self._eng, self._children)

    def as_dict(self):
        return {
                key: value
                for key, value in zip(self._fields, self.as_tuple())
                }

    def set_children(self, children=None):
        if children:
            #it is a tuple can be only changed by append
            self._children.extend(children)

    def __int__(self):
        """
            int(Parameter) will return the raw value of the parameter
        """
        return self.get_raw_int()

    def __str__(self):
        return str(self.parameter())
    #return self.get(key)

    def __setitem__(self, key, value):
        if key == 'name':
            self._name = value
        elif key == 'raw':
            self._raw = value
        elif key == 'eng':
            self._eng = value
        elif key == 'children':
            self._children = value
        elif key in ['error','aux']:
            self._aux= value
        else:
            raise KeyError(key)

    def __getattr__(self, key):
        """
         support parameter.key to access the value
        """
        if key == 'name':
            return self._name
        elif key == 'raw':
            return self._raw
        elif key == 'eng':
            return self._eng
        elif key == 'children':
            return self._children
        elif key == 'raw_int':
            return self.get_raw_int()
        elif key == 'desc':
            return STIX_IDB.get_parameter_description(self._name)
        elif key == 'param':
            return self.as_tuple()
        elif key in ['error','aux']:
            return self._aux
        else:
            raise KeyError(key)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.__getattr__(key)
        elif isinstance(key, int):
            return self.as_tuple[key]
        else:
            raise KeyError(key)


class Packet(object):
    def __init__(self, a=None, b=None, deep_copy=False):
        self._header = {}
        self._parameters = []
        self._current_node_name = '/'

        if a is None:
            a = {}
        if b is None:
            b = {}
        if 'header' in a and 'parameters' in a:
            self._header = copy_object(a['header'], deep_copy)
            self._parameters = copy_object(a['parameters'], deep_copy)
        elif isinstance(a, dict):
            self._header = copy_object(a, deep_copy)
            if isinstance(b, list):
                self._parameters = copy_object(b, deep_copy)

    def is_valid(self):
        if self._header and self._parameters:
            return True
        return False

    def get_raw_length(self):
        if not self._header:
            return 0
        TMTC = self._header['TMTC']
        length_header = self._header['length']
        if TMTC == 'TM':
            return length_header - 9 + 16
        return length_header + 1 - 4 + 10

    def as_dict(self):
        return {'header': self._header, 'parameters': self._parameters}

    def __str__(self):
        return pprint.pformat(self.as_dict())

    def __getattr__(self, attr):
        if attr == 'header':
            return self._header
        elif attr == 'parameters':
            return self._parameters
        elif attr in self._header:
            return self._header['{}'.format(attr)]
        elif attr == 'raw_len':
            return self.get_raw_length()
        else:
            return self.get_nodes(attr)

    def walk_parameters(self, func):
        """
        walk parameters, and apply function to it
        """
        pass


    def get_nodes(self, pattern, parameters=None):
        if pattern.startswith('/'):
            pattern = pattern[1:]
        fields = pattern.split('/')

        results = []
        if not fields:
            return []
        if not parameters:
            parameters = self._parameters
        try:
            field = fields.pop(0)
        except IndexError:
            return []
        for e in parameters:
            param = Parameter(e)
            if param['name'] in field or '*' in field:
                if fields:
                    ret = self.get_nodes('/'.join(fields), param['children'])
                    if ret:
                        results.append(ret)
                else:

                    results.append(param.as_tuple())

        return results

    def test_conditions(self, parameter, conditions):

        if not conditions:
            return True
        eval_str = conditions
        nmax = 20
        if 'raw' in conditions:
            eval_str = eval_str.replace('raw', str(parameter['raw']), nmax)
        if 'eng' in conditions:
            eval_str = eval_str.replace('eng', str(parameter['eng']), nmax)

        try:
            return eval(eval_str)
        except Exception as e:
            logger.error("test condition:" + str(e))
            return True

    def get(self, pattern, conditions='', parameters=None):
        """
          pattern examples:
              pattern='NIX00159/NIX00146.eng'
                  return the eng. values of all NIX00146 under NIX00159
              pattern='NIX00159/NIX00146.raw/*.eng'
                  return the children's eng. value of all NIX00146
            examples:
            p2.get('NIXG0020/*.name','raw <0 and eng >0')
            p2.get('NIXG0020/*.raw','raw <0 and eng >0')

        """
        if pattern.startswith('/'):
            pattern = pattern[1:]

        fields = pattern.split('/')
        results = []
        if not fields:
            return []
        if not parameters:
            parameters = self._parameters
        try:
            field = fields.pop(0)
        except IndexError:
            return []
        for e in parameters:
            param = Parameter(e)
            if param['name'] in field or '*' in field:
                if fields:
                    ret = self.get('/'.join(fields), conditions,
                                   param['children'])
                    if ret:
                        results.append(ret)
                else:
                    if not self.test_conditions(param, conditions):
                        continue

                    if field.endswith('eng'):
                        results.append(param['eng'])
                    elif field.endswith('name'):
                        results.append(param['name'])
                    else:
                        results.append(param['raw_int'])
        return results

    def cd(self, node_name=None):

        if not node_name or node_name == '/' or node_name == '':
            self._current_node_name = '/'

        else:
            self._current_node_name = os.path.abspath(
                os.path.join(self._current_node_name, node_name))

        print('current path: {}'.format(self._current_node_name))

    def __getitem__(self, key):
        #access parameter like packet['parameter']
        if isinstance(key, str):
            return self.__getattr__(key)
        elif isinstance(key, int):
            return Parameter(self._parameters[key])

    def isa(self, SPIDs):

        if self._header.get('TMTC', None) == 'TC':
            return False
        if not isinstance(SPIDs, list):
            SPIDs = [int(SPIDs)]
        try:
            if int(self._header['SPID']) in SPIDs:
                return True
        except (KeyError, ValueError) as e:
            logger.warn("SPID not found in the packet")
        return False

    def ls(self, node_name=None):
        if not node_name:
            node_name = self._current_node_name
        num_fields = len(node_name.split('/'))
        pprint.pprint(self.get_nodes(node_name))

    @staticmethod
    def merge(packets, SPIDs, value_type='raw', remove_duplicates=True):
        if not isinstance(SPIDs, list):
            SPIDs = [SPIDs]

        result = {}
        hash_list = []
        min_id = math.inf
        max_id = -math.inf
        pkt_len = 0
        
        seg_flag_counter={0:0,1:0,2:0,3:0}
        
        for pkt in packets:
            try:
                if int(pkt['header']['SPID']) not in SPIDs:
                    continue
            except (ValueError, KeyError):
                continue

            if remove_duplicates:
                if pkt['hash'] in hash_list:
                    logger.info(f'Duplicated packet {pkt["_id"]} excluded!')
                    continue
                hash_list.append(pkt['hash'])

            min_id, max_id = min(min_id, pkt['_id']), max(max_id, pkt['_id'])
            pkt_len += 1

            if set(SPIDs).intersection([54115, 54143]): 
                #l1 or L4 report, Hualin 2022-05-10, fix the issues of the same data were requested twice (included into different IORs)
                seg_flag_counter[pkt['header']['seg_flag']] += 1
                if seg_flag_counter[1] >1 or seg_flag_counter[2]>1:
                    logger.warning(f'{pkt["_id"]} ignored, might be a duplicated requests')
                    continue
            Packet.merge_headers(result, pkt['header'])
            Packet.merge_parameters(result, pkt['parameters'], value_type)

        result['min_id'] = min_id
        result['max_id'] = max_id
        result['num_packets'] = pkt_len

        return result

    @staticmethod
    def merge_headers(result, header):
        """
        merge header
        """
        for key, value in header.items():
            if key not in result:
                result[key] = [value]
            else:
                result[key].append(value)

    @staticmethod
    def merge_parameters(result, parameters, value_type):
        """merge packets recursively 
          
          history:
          2020-05-23, removed conversion of parameter to Parameter object

        """

        if not parameters:
            return
        for p in parameters:
            #param = Parameter(p)
            name, raw, eng, children=p

            value = eng if (value_type == 'eng' and  eng) else raw

            #if 'NIXG' not in name:

            if name not in result:
                result[name] = []
            result[name].append(value)
            if children:
                Packet.merge_parameters(result, children, value_type)
            #doesn't work

    def children_as_dict(self, parameter_names=None, children=None):
        if not children:
            children = self._parameters
        param_dict = {}
        for e in children:
            param = Parameter(e)
            if 'NIXG' in param['name']:
                continue
            if parameter_names:
                if param['name'] not in parameter_names:
                    continue
            if param['name'] in param_dict:
                param_dict[param.name].append(param)
            else:
                param_dict[param.name] = [param]

        return param_dict

    def index(self, parameter_name):
        #get parameter index
        #only looks for parameters whose depth == 0
        for i, e in enumerate(self._parameters):
            if e[0] == parameter_name:
                return i
        return -1

    def get_one(self, parameter_name, parameters=None):
        #get the first parameter
        if parameters == None:
            parameters = self._parameters
        for e in parameters:
            if e[0] == parameter_name:
                return e
            if e[3]:
                ret = self.get_one(parameter_name, e[3])
                if ret:
                    return ret

        return None

    def get_many(self, selectors, parent=[], parameters=None):
        """
        Extract parameter values from packet
          supports multiple selectors 

          selectors = ['NIX00159/NIX00146.eng', 'NIXG0020/*.name']
          pattern examples:
              pattern='NIX00159/NIX00146.eng'
                  return the eng. values of all NIX00146 under NIX00159
              pattern='NIX00159/NIX00146.raw/*.eng'
                  return the children's eng. value of all NIX00146
            examples:
            p2.get('NIXG0020/*.name','raw <0 and eng >0')
            p2.get('NIXG0020/*.raw','raw <0 and eng >0')

        """
        if not parameters:
            parameters = self._parameters

        for e in parameters:
            param = Parameter(e)
            parent.append(param['name'])

            if param['children'] and any_list_contains(selectors, parent):
                self.get_many(selectors, parent, param['children'])
        '''
            
            if param['name'] in field or '*' in field:
                if fields:
                    ret = self.get('/'.join(fields), conditions,
                                   param['children'])
                    if ret:
                        results.append(ret)
                else:
                    if not self.test_conditions(param, conditions):
                        continue

                    if field.endswith('eng'):
                        results.append(param['eng'])
                    elif field.endswith('name'):
                        results.append(param['name'])
                    else:
                        results.append(param['raw_int'])
        return results
        '''
