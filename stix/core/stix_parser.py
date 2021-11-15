#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : parser.py
# @date         : Feb. 11, 2019
# @author       : Hualin Xiao (hualin.xiao@fhnw.ch)
# @description:
#               STIX telemetry raw data parser
# @TODO
#              add checks for CRC for TC packets
import os
import math
import re
import struct as st
import binascii
import hashlib
import pathlib
import xmltodict
from scipy import interpolate
from dateutil import parser as dtparser
from stix.core import stix_header
from stix.core import stix_idb
from stix.core import stix_global
from stix.core import stix_logger
from stix.core import stix_writer
from stix.core import stix_context
from stix.core.stix_datatypes import Parameter
from stix.core import stix_decompressor
from stix.spice import stix_datetime

CONTEXT_UNPACK_FORMAT = ['B', '>H', 'BBB', '>I']
UNSIGNED_UNPACK_FORMAT = ['B', '>H', 'BBB', '>I', 'BBBBB', '>IH']
SIGNED_UNPACK_FORMAT = ['b', '>h', 'bbb', '>i', 'bbbbb', '>ih']
HEX_SPACE = '0123456789ABCDEFabcdef'
SCET_PARAMETERS = [
    'NIX00402', 'NIX00445', 'NIX00287', 'PIX00455', 'PIX00456', 'PIX0021',
    'PIX0022', 'PIX0025', 'PIX0026', 'PIX00086', 'PIX00087', 'PIX00009'
]

PARMETERS_CALIBRATION_ENABLED = [
    'NIX00101', 'NIX00102', 'NIX00125', 'NIXD0003'
]
PARMETERS_CALIBRATION_ENABLED.extend(SCET_PARAMETERS)

STIX_IDB = stix_idb.stix_idb()
STIX_DECOMPRESSOR = stix_decompressor.StixDecompressor()
logger = stix_logger.get_logger()


def detect_filetype(filename):
    """detect file type when its type is not specified
    Parameters:
       filename: input filename
    Returns:
       file type as a string
    """
    filetype = None
    extension = pathlib.Path(filename).suffix
    ext = extension[1:]

    if ext in ('ascii', 'bin', 'hex'):
        return ext
    elif ext in ('binary', 'raw', 'BDF', 'dat'):
        return 'bin'
    try:
        fin = open(filename, 'r')
        buf = fin.read(1024).strip()
        if ext == 'xml':
            if 'PktTcReportResponse' in buf:
                filetype = 'xml.tc.hist'
            elif 'PktRawResponse' in buf:
                filetype = 'xml.tm.raw'
        if not filetype:
            data = re.sub(r"\s+", "", buf)
            filetype = 'hex'
            for c in data:
                if c not in HEX_SPACE:
                    filetype = 'ascii'
                    #possibly ASCII
                    break
    except UnicodeDecodeError:
        filetype = 'bin'
    finally:
        fin.close()
    return filetype


def slice_bits(data, offset, num_bits):
    """slice bits from buffer 
    Parameters:
        data: input integer
        offset: offset in units of bits
        num_bits:  number of bits to extract
    """
    return (data >> offset) & ((1 << num_bits) - 1)


def to_int(value):
    try:
        return int(value)
    except (TypeError, IndexError, ValueError):
        return None


def unpack_integer(raw, structure):
    """unpack variables from an integer 

    Parameters:
        raw: input integer 
        structure: variable structure

    Returns:
      a dictionary containing extracted variables
    """
    result = {}
    for name, bits in structure.items():
        result[name] = slice_bits(raw, bits[0], bits[1])
    return result


def get_from_bytearray(buf, i, width=1):
    """slice bytes from a byte array

    """
    data = buf[i:i + width]
    length = len(data)
    if length != width:
        return False, i + length, data
    return True, i + length, data


def find_next_header(buf, i):
    """find next TC/TM packet header
    Parameters:
        buf: buffer
        i: start position
    Returns:
        index or EOF
    """
    length = len(buf)
    while i < length:
        x = buf[i]
        if x in stix_header.HEADER_FIRST_BYTE:
            return i
        else:
            i += 1
    return stix_global.EOF


class StixParameterParser(object):
    def __init__(self):
        pass

    def decode_buffer(self,
                      in_data,
                      param_type,
                      offset,
                      offset_bits,
                      length,
                      param_name=''):
        """Decode a parameter  from  buffer
        Parameters:
            parameter_type: parameter type
            offset: offset in units of bits
            offset_bits: offset within a byte
        Returns:
             a tuple with decoded parameter value 

        """
        nbytes = math.ceil((length + offset_bits) / 8.)

        raw_bin = in_data[int(offset):int(offset + nbytes)]

        if nbytes != len(raw_bin):
            logger.error(
                'Parameter {} length mismatch.  Expect: {} real: {}'.format(
                    param_name, nbytes, len(raw_bin)))
            return ''
        bin_struct = str(nbytes) + 's'  #signed char
        if param_type == 'U' and nbytes <= 6:
            bin_struct = UNSIGNED_UNPACK_FORMAT[nbytes - 1]
        elif param_type == 'I' and nbytes <= 6:
            bin_struct = SIGNED_UNPACK_FORMAT[nbytes - 1]
        elif param_type == 'T':
            bin_struct = '>IH'
        elif param_type == 'CONTEXT' and nbytes <= 4:
            bin_struct = CONTEXT_UNPACK_FORMAT[nbytes - 1]

        raw = st.unpack(bin_struct, raw_bin)

        if len(raw) == 1:
            result = raw[0]
            if length == 8 and offset_bits == 0:
                #only used to extract archive memory information
                return result
            elif length < 16:  #and length % 8 != 0:
                # bit-offset only for 8 bits or 16 bits integer
                # only valid for STIX
                start_bit = nbytes * 8 - (offset_bits + length)
                return slice_bits(raw[0], start_bit, length)
            return result
        elif len(raw) == 2:
            if param_type == 'T':
                return round(float(raw[0]) + float(raw[1]) / 65536., 3)
            logger.warning(
                'Invalid unpacking parameter type: {}'.format(param_type))
            return raw_bin
        elif len(raw) == 3:  # 24-bit integer, a timestamp probably
            value = (raw[0] << 16) | (raw[1] << 8) | raw[2]
            if length < 16 and length % 8 != 0:
                start_bit = nbytes * 8 - (offset_bits + length)
                value = slice_bits(value, start_bit, length)
            return value

        return raw_bin

    def raw_to_eng(self, param_name, ref, raw_value, tmtc='TM'):
        """ convert raw value to engineering value

        Parameters:
            param_name: parameter name
            ref:  parameter reference name as defined in IDB
            parameter_type: parameter type as defined in IDB
            raw: raw value
            tmtc: string 'TM' or 'TC'

        Return:
           engineering value or an empty string if its has no engineering value
        """

        if raw_value is None:
            return ''
        if param_name in SCET_PARAMETERS:
            #convert SCET to UTC
            return stix_datetime.scet2utc(int(raw_value))

        if tmtc == 'TC':
            if ref:
                return STIX_IDB.tcparam_interpret(ref, raw_value)
            return ''
        if param_name == 'NIX00125':  #temperature calibration factors, parameter from Richard on June 22, 2020
            try:
                x = math.log(12000. * raw_value / (4095 - raw_value))
                return round(
                    1 /
                    (9.27e-4 + 2.34e-4 * x + 9.09e-7 * x**2 + 1.04e-7 * x**3) -
                    273.15, 2)
            except ValueError:
                logger.warning(
                    'Could not calibrate NIX00125 temperature raw value :{}'.
                    format(raw_value))
        if param_name == 'NIXD0003':
            return round(raw_value / 2.5, 1)  # archive memory to raw value

        #conversion based on the equation in SIRIUS source code
        #    return (raw_value * 1.1 * 3.0 / 4095 - 1.281) * 213.17
        #elif param_name == 'NIX00102':
        #temperature std. deviations
        #    return (raw_value * 1.1 * 3.0 / 4095) * 213.17

        if not ref:
            return ''
        prefix = re.split(r'\d+', ref)[0]

        if prefix in ['CIXTS', 'CAAT', 'CIXT']:
            # textual interpret
            rows = STIX_IDB.textual_interpret(ref, raw_value)
            if rows:
                return rows[0][0]
            logger.warning(
                'Missing textual calibration info. for {}'.format(ref))
            return ''
        elif prefix == 'CIXP':
            # calibration
            # Ref SCOS-2000 Database Import ICD
            rows = STIX_IDB.get_calibration_curve(ref)
            if len(rows) <= 1:
                logger.warning(
                    'Invalid calibration parameter {}: at least two data points needed '
                    .format(ref))
                return ''
            elif len(rows) >= 2:
                x_points = [float(row[0]) for row in rows]
                y_points = [float(row[1]) for row in rows]

                if len(rows) == 2:
                    return round((y_points[1] - y_points[0]) /
                                 (x_points[1] - x_points[0]) *
                                 (raw_value - x_points[0]) + y_points[0], 3)
                try:
                    tck = interpolate.splrep(x_points, y_points)
                    val = interpolate.splev(raw_value, tck)
                    ret = round(float(val), 3)
                except Exception as e:
                    logger.warning('Failed to calibrate {} due to {}'.format(
                        ref, str(e)))
                    ret = ''
                return ret
            logger.warning('Missing calibration factors for {}'.format(ref))
            return ''
        elif prefix == 'NIX':
            logger.warning('No information to convert {} . '.format(ref))
            return ''
        elif prefix == 'CIX':
            rows = STIX_IDB.get_calibration_polynomial(ref)
            if rows:
                pol_coeff = ([float(x) for x in rows[0]])
                x_points = ([math.pow(raw_value, i) for i in range(0, 5)])
                sum_value = 0
                for coeff, xval in zip(pol_coeff, x_points):
                    sum_value += coeff * xval
                return round(sum_value, 3)
            logger.warning('Missing calibration factors for {}'.format(ref))
            return ''

        return ''

    def decode_parameter(self,
                         buf,
                         name,
                         offset,
                         offset_bits,
                         width,
                         ptc,
                         pfc,
                         cal_ref='',
                         tmtc='TM',
                         calibration_enabled=True):
        """Decode and calibrate a parameter 

        Parameters:
            buf:  buffer
            name: parameter name 
            offset: offset relative to the beginning of the buffer
            offset_bit: offset relative to the beginning of the integer
            width: width in units of bits
            ptc:  parameter type code as defined   and SCOS 2000
            pfc:  parameter format code as defined and SCOS 2000
            cal_ref: calibration reference name
            tmtc:  packet type 'TM' or 'TC'
            calibration_enabled: calibrate the parameter or not
        Returns:
           parameter in STIX  parameter type. 

        """
        param_type = STIX_IDB.get_s2k_parameter_types(ptc, pfc)  #1.79
        raw_value = self.decode_buffer(buf, param_type, offset, offset_bits,
                                       width, name)
        eng_value = ''
        if calibration_enabled:
            eng_value = self.raw_to_eng(name, cal_ref, raw_value, tmtc)

        if raw_value != '':
            #try to decompress raw value
            raw_int = raw_value
            try:
                raw_int = int(raw_int)
                result = STIX_DECOMPRESSOR.decompress_raw(name, raw_int)
                if result is not None:
                    eng_value = result
            except (TypeError, ValueError):
                pass

        #param = Parameter(
        return (name, raw_value, eng_value, [])
        #return param['param']


class StixVariableTelemetryPacketParser(StixParameterParser):
    """Variable length telemetry packet parser
    """
    def __init__(self):
        super(StixVariableTelemetryPacketParser, self).__init__()
        self.last_spid = -1
        self.last_num_bits = 0
        self.last_data_width = 0
        self.current_offset_bit = 0
        self.length_min = 0
        self.nodes = []
        self.buffer = None
        self.spid = 0
        self.current_offset = 0
        self.last_offset = 0
        self.results_tree = []
        self.results_dict = {}

    def init_nodes(self):
        """Initialize parse tree
        """
        self.nodes = []
        self.nodes.append(self.create_parse_node('top', counter=1))
        self.nodes[0]['children'] = []
        self.length_min = 0

    def parse(self, data, spid):
        """parse a variable TM packet
        """

        self.last_data_width = 0
        self.last_num_bits = 0
        #width of last parameter

        self.buffer = data
        self.spid = spid

        self.current_offset = 0
        self.last_offset = 0
        self.current_offset_bit = 0
        self.results_tree = []
        self.results_tree[:] = []
        self.results_dict = {}
        self.results_dict.clear()
        if spid != self.last_spid:
            self.init_nodes()
            self.build_parse_tree()

        self.last_spid = spid
        packet_length = len(self.buffer)
        if self.length_min > packet_length:
            return 0, None, stix_global.VARIABLE_PACKET_LENGTH_MISMATCH
        self.walk(self.nodes[0], self.results_tree)
        return self.current_offset, self.results_tree, stix_global.OK

    def create_parse_node(self,
                          name,
                          parameter=None,
                          counter=0,
                          children=None):

        if children is None:
            children = []
        node = {
            'name': name,
            'counter': counter,
            'parameter': parameter,
            'children': children
        }
        return node

    def register_parameter(self, mother, param_PCF):
        children = []
        name = param_PCF['PCF_NAME']
        node = self.create_parse_node(name, param_PCF, 0, children)
        width = param_PCF['PCF_WIDTH']
        if width % 8 == 0:
            self.length_min += width / 8
        mother['children'].append(node)
        return node

    def build_parse_tree(self):
        """
        To build a parameter parse tree
        """
        param_pcf_structures = STIX_IDB.get_variable_packet_structure(
            self.spid)

        mother = self.nodes[0]
        repeater = [{
            'node': mother,
            'counter': stix_global.MAX_NUM_PARAMETERS
        }]
        for par in param_pcf_structures:
            if repeater:
                for e in reversed(repeater):
                    e['counter'] -= 1
                    if e['counter'] < 0:
                        repeater.pop()
                        #root will be never popped
            mother = repeater[-1]['node']
            node = self.register_parameter(mother, par)

            rpsize = par['VPD_GRPSIZE']
            if rpsize > 0:
                mother = node
                repeater.append({'node': node, 'counter': rpsize})

    def walk(self, mother, parameters):
        if not mother:
            return
        counter = mother['counter']
        #mother_name = mother['name']
        for i in range(0, counter):
            for pnode in mother['children']:
                if not pnode or self.current_offset > len(self.buffer):
                    return
                ret = self.parse_node(pnode)
                #param = Parameter(ret)
                param = ret
                if pnode['children']:
                    #num_children = param.raw_int
                    num_children = to_int(param[1])
                    is_valid = False
                    if isinstance(num_children, int):
                        if num_children > 0:
                            pnode['counter'] = num_children
                            is_valid = True
                            self.walk(pnode, param[3])
                            # if mother_name == 'NIXD0159':
                            # correction for Science L0 data, ICD not fully consistent with IDB,
                            # STIX ICD-0812-ESC Table 93 P123
                            # self.correct_science_L0(param)

                    if not is_valid:
                        if param[0] != 'NIXD0159':
                            #repeater NIXD0159 can be zero according to STIX ICD-0812-ESC Table 93 P123
                            logger.warning(
                                'Repeater {}  has an invalid value: {}'.format(
                                    pnode['name'], param[1]))

                parameters.append(param)

    def parse_node(self, node):
        """
        decode a parameter.
        The offset, offset bit and parameter name are described in 'node'
        """
        par = node['parameter']

        offset_bits = int(par['VPD_OFFSET'])
        width = int(par['PCF_WIDTH'])
        ptc = int(par['PCF_PTC'])
        pfc = int(par['PCF_PFC'])
        cal_ref = par['PCF_CURTX']
        name = node['name']

        if width % 8 != 0:
            if offset_bits < 0:
                self.current_offset_bit = self.last_data_width + offset_bits
            else:
                self.current_offset_bit += self.last_num_bits + offset_bits
            self.last_num_bits = width
        elif width % 8 == 0:
            self.current_offset_bit = 0
            self.last_offset = self.current_offset
            self.current_offset += width / 8
            self.last_num_bits = 0
            self.last_data_width = width

        offset = self.last_offset
        offset_bits = self.current_offset_bit
        calibration_enabled = False
        if name in PARMETERS_CALIBRATION_ENABLED:
            calibration_enabled = True

        #Don't convert raw to eng for variable length packets except those in the list
        return self.decode_parameter(self.buffer, name, offset, offset_bits,
                                     width, ptc, pfc, cal_ref, 'TM',
                                     calibration_enabled)


class StixContextParser(StixParameterParser):
    '''Context file parser
       Context file structure is not described in STIX IDB
    '''
    def __init__(self):
        super(StixContextParser, self).__init__()

    def parse(self, buf):
        #parsing based on the FSW source code  ContextMgmt

        offset = 0
        #offset in units of bits
        parameters = []
        param_id = 0
        for name, width in stix_context.CONTEXT_PARAMETER_BIT_SIZE:
            #width also in units of bits
            offset_bytes = int(offset / 8)
            offset_bits = offset % 8
            children = []
            raw_values = None
            if name in stix_context.ASIC_REGISTERS:
                children = self.parse_asic_registers(buf, offset)
                raw_values = (len(children), )  #as a repeater
            else:
                raw_values = self.decode_buffer(buf, 'CONTEXT', offset_bytes,
                                                offset_bits, width)
            parameters.append((name, raw_values, '', children))

            offset += width

            param_id += 1
        return parameters

    def parse_asic_registers(self, buf, offset):
        parameters = []
        for name, width in stix_context.CONTEXT_REGISTER_BIT_SIZE:
            offset_bytes = int(offset / 8)
            offset_bits = offset % 8
            raw_value = self.decode_buffer(buf, 'CONTEXT', offset_bytes,
                                           offset_bits, width)
            offset += width
            if raw_value is not None and raw_value != '':
                param = (stix_context.CONTEXT_REGISTER_DESC[name], raw_value,
                         '', [])
                parameters.append(param)
        return parameters


class StixTelecommandParser(StixParameterParser):
    """
        STIX telecommand packet parser
    """
    def __init__(self):
        super(StixTelecommandParser, self).__init__()
        self.length_min = 0
        self.nodes = []
        self.last_tc_name = ''
        self.tc_name = ''
        self.current_bit_offset = 0
        self.results_tree = []
        self.param_structure = []
        self.buffer = None

    def parse(self, name, buf):
        """
            To parse a telecommand
        """
        self.buffer = buf
        self.tc_name = name
        self.param_structure = []
        self.current_bit_offset = 0
        is_variable = STIX_IDB.is_variable_length_telecommand(name)
        self.param_structure = STIX_IDB.get_telecommand_structure(name)
        if is_variable:
            return self.parse_variable_telecommand()
        return self.parse_fixed_telecommand()

    def parse_fixed_telecommand(self):
        params = []
        for par in self.param_structure:
            result = self.parse_one(par,
                                    is_fixed=True,
                                    calibration_enabled=True)
            params.append(result)
        return int(self.current_bit_offset / 8), params, stix_global.OK

    def parse_one(self, par, is_fixed=True, calibration_enabled=True):
        cdf_type = par['CDF_ELTYPE']
        param_name = ''
        ptc = 0
        pfc = 0
        if cdf_type == 'A':
            #fixed area SCOS-2000 ICD page 63
            param_name = par['CDF_DESCR']
            ptc = 3  # works for STIX
            pfc = 4
        else:
            ptc = int(par['CPC_PTC'])
            pfc = int(par['CPC_PFC'])
            param_name = par['CDF_PNAME']
        cal_ref = par['CPC_PAFREF']
        width = int(par['CDF_ELLEN'])

        if is_fixed:
            offset = int(int(par['CDF_BIT']) / 8)
            offset_bits = int(par['CDF_BIT']) % 8
            self.current_bit_offset = int(par['CDF_BIT']) + width
        else:
            #CDF_ELLIEN is not respected  for variable telecommands
            offset = int(self.current_bit_offset / 8)
            offset_bits = int(self.current_bit_offset % 8)
            # need to be checked, may be wrong
            self.current_bit_offset += width

        parameter = self.decode_parameter(
            self.buffer,
            param_name,
            offset,
            offset_bits,
            width,
            ptc,
            pfc,
            cal_ref,
            'TC',
            calibration_enabled=calibration_enabled)
        return parameter

    def parse_variable_telecommand(self):
        self.current_bit_offset = 0
        self.length_min = 0
        self.results_tree = []
        self.nodes = []
        self.nodes.append(self.create_parse_node('top', counter=1))
        self.nodes[0]['children'] = []
        self.length_min = 0
        self.build_parse_tree()
        packet_length = len(self.buffer)
        if self.length_min > packet_length:
            return 0, None, stix_global.VARIABLE_PACKET_LENGTH_MISMATCH
        self.walk(self.nodes[0], self.results_tree)
        return int(self.current_bit_offset /
                   8), self.results_tree, stix_global.OK

    def create_parse_node(self,
                          name,
                          parameter=None,
                          counter=0,
                          children=None):
        if children is None:
            children = []
        node = {
            'name': name,
            'counter': counter,
            'parameter': parameter,
            'children': children
        }
        return node

    def register_parameter(self, mother, idb_param):
        children = []
        name = idb_param['CDF_PNAME']
        node = self.create_parse_node(name, idb_param, 0, children)
        width = idb_param['CDF_ELLEN']
        if width % 8 == 0:
            self.length_min += width / 8
        mother['children'].append(node)
        return node

    def build_parse_tree(self):
        """
        To build a parameter parse tree
        """
        mother = self.nodes[0]
        repeater = [{
            'node': mother,
            'counter': stix_global.MAX_NUM_PARAMETERS
        }]
        for par in self.param_structure:
            if repeater:
                for e in reversed(repeater):
                    e['counter'] -= 1
                    if e['counter'] < 0:
                        repeater.pop()
            mother = repeater[-1]['node']
            node = self.register_parameter(mother, par)
            rpsize = par['CDF_GRPSIZE']
            if rpsize > 0:
                mother = node
                repeater.append({'node': node, 'counter': rpsize})

    def walk(self, mother, parameters):
        if not mother:
            return
        counter = mother['counter']
        for i in range(0, counter):
            for pnode in mother['children']:
                if not pnode or self.current_bit_offset > 8 * len(self.buffer):
                    return
                ret = self.parse_one(pnode['parameter'],
                                     False,
                                     calibration_enabled=True)
                param = ret
                if pnode['children']:
                    #num_children = param['raw_int']
                    num_children = to_int(param[1])
                    is_valid = False
                    if isinstance(num_children, int):
                        if num_children > 0:
                            pnode['counter'] = num_children
                            is_valid = True
                            #self.walk(pnode, param.children)
                            self.walk(pnode, param[3])
                    if not is_valid:
                        logger.warning(
                            'Repeater {}  has an invalid value: {}. '.format(
                                pnode['name'], param[1]))

                #parameters.append(param.as_tuple())
                parameters.append(param)


class StixTCTMParser(StixParameterParser):
    def __init__(self):
        super(StixTCTMParser, self).__init__()
        self.vp_tm_parser = StixVariableTelemetryPacketParser()
        self.tc_parser = StixTelecommandParser()
        self.context_parser = StixContextParser()
        self.selected_services = []
        self.selected_spids = []
        #self.store_binary = True
        #self.is_live_hex_stream = False
        #self.decoded_packets = []
        self.raw_filename = ''
        self.in_filesize = 0
        self.stop_parsing = False

        self.receipt_utc = ''

        self.store_packet_enabled = True
        self.packet_writer = None
        self.S20_excluded = False

        self.stix_alerts = []
        self.parser_counter = {}
        #parser counter
        self.reset_counter()

    def reset_counter(self):
        self.parser_counter = {
            'num_tm': 0,
            'num_tc': 0,
            'num_tm_parsed': 0,
            'num_tc_parsed': 0,
            'num_bad_bytes': 0,
            'num_bad_headers': 0,
            'total_length': 0,
            'num_filtered': 0,
            'spid': {}
        }

    def inc_counter(self, t, value=1):
        #increase count t by the given value or increase spid counters
        if t == 'spid':
            if str(value) not in self.parser_counter['spid']:
                self.parser_counter['spid'][str(value)] = 0
            self.parser_counter['spid'][str(value)] += 1
            #mongodb only supports string as keys
            return
        self.parser_counter[t] += value

    def get_summary(self):
        return self.parser_counter

    def get_stix_alerts(self):
        # return a summary of TM(5, 2..4)  and TM(1,2), TM(1,8)
        return self.stix_alerts

    def kill(self):
        self.stop_parsing = True

    def exclude_S20(self):
        self.S20_excluded = True

    def set_packet_buffer_enabled(self, status):
        """
            store packets in a list
        """
        self.store_packet_enabled = status

    def set_packet_filter(self, selected_services=None, selected_spids=None):
        """ only decoded packets with the given services or spids
        """
        if selected_spids is None:
            selected_spids = []
        if selected_services is None:
            selected_services = []

        self.selected_services = selected_services
        self.selected_spids = selected_spids

    #def set_store_binary_enabled(self, status):
    #    """
    #      store raw binary  in the output
    #    """
    #    self.store_binary = status

    def parse_telemetry_header(self, packet):
        """ see STIX ICD-0812-ESC  (Page # 57) """
        if len(packet) < 16:
            logger.warning("Packet length < 16. Packet ignored!")
            return stix_global.PACKET_TOO_SHORT, None
        if ord(packet[0:1]) not in stix_header.TM_HEADER_FIRST_BYTE:
            return stix_global.HEADER_FIRST_BYTE_INVALID, None
        header_raw = st.unpack('>HHHBBBBIH', packet[0:16])
        header = {}
        for h, s in zip(header_raw, stix_header.TELEMETRY_RAW_STRUCTURE):
            header.update(unpack_integer(h, s))
        status = self.check_header(header, 'tm')
        if status == stix_global.OK:
            header['segmentation'] = stix_header.PACKET_SEG[header['seg_flag']]
            header['TMTC'] = 'TM'
            header[
                'SCET'] = header['fine_time'] / 65536. + header['coarse_time']
            header['raw_length'] = header['length'] - 9 + 16
        return status, header

    def check_header(self, header, tmtc='tm'):
        # header validation
        constrains = None
        if tmtc == 'tm':
            constrains = stix_header.TELEMETRY_HEADER_CONSTRAINTS
        else:
            constrains = stix_header.TELECOMMAND_HEADER_CONSTRAINTS
        for name, lim in constrains.items():
            if header[name] not in lim:
                logger.warning(
                    'Header {} value {} violates the range: {} '.format(
                        name, header[name], lim))
                return stix_global.HEADER_INVALID
        return stix_global.OK

    def parse_data_field_header(self, header, data, length):
        """ Decode the data field header
        """
        service_type = header['service_type']
        service_subtype = header['service_subtype']
        offset, width = STIX_IDB.get_packet_type_offset(
            service_type, service_subtype)
        # see solar orbit ICD Page 36
        ssid = -1
        if offset != -1:
            start = offset - 16  # 16bytes read already
            end = start + width / 8  # it can be : 0, 16,8
            bin_struct = '>B'
            if width == 16:
                bin_struct = '>H'
            raw = data[int(start):int(end)]
            if not raw:
                return stix_global.PACKET_TOO_SHORT
            res = st.unpack(bin_struct, raw)
            ssid = res[0]
        info = STIX_IDB.get_packet_type_info(service_type, service_subtype,
                                             ssid)
        if not info:
            return stix_global.NO_PID_INFO_IN_IDB
        header['descr'] = info['PID_DESCR']
        header['SPID'] = info['PID_SPID']
        header['TPSD'] = info['PID_TPSD']
        header['length'] = length
        header['SSID'] = ssid
        return stix_global.OK

    def parse_fixed_telemetry_packet(self, buf, spid):
        """ Extract parameters for a fixed data packet 
            see Solar orbit IDB.ICD section 3.3.2.5.1
        """
        if spid == 54331:
            #context file report parsing
            # not to use IDB
            return self.context_parser.parse(buf)
        parameters = []
        param_structures = STIX_IDB.get_fixed_packet_structure(spid)
        for par in param_structures:
            offset = int(par['PLF_OFFBY']) - 16
            offset_bits = int(par['PLF_OFFBI'])
            width = int(par['PCF_WIDTH'])
            ptc = int(par['PCF_PTC'])
            pfc = int(par['PCF_PFC'])
            name = par['PCF_NAME']
            cal_ref = par['PCF_CURTX']
            param = self.decode_parameter(buf, name, offset, offset_bits,
                                          width, ptc, pfc, cal_ref, 'TM', True)
            parameters.append(param)
        return parameters

    def parse_telecommand_header(self, buf, ipos):
        # see STIX ICD-0812-ESC  (Page 56)
        if buf[ipos] not in stix_header.TC_HEADER_FIRST_BYTE:
            return stix_global.HEADER_FIRST_BYTE_INVALID, None
        try:
            header_raw = st.unpack('>HHHBBBB', buf[ipos:ipos + 10])
        except Exception as e:
            logger.error(str(e))
            return stix_global.HEADER_RAW_LENGTH_VALID, None, None
        header = {}
        for h, s in zip(header_raw, stix_header.TELECOMMAND_RAW_STRUCTURE):
            header.update(unpack_integer(h, s))
        status = self.check_header(header, 'tc')

        header['raw_length'] = header['length'] - 3 + 10

        subtype = -1
        if (header['service_type'], header['service_subtype']) in [(237, 7),
                                                                   (236, 6)]:
            #extra bits are needed to identify TC type, see ICD
            try:
                subtype = st.unpack('B', buf[ipos + 10:ipos + 11])[0]
                header['subtype'] = subtype
            except Exception as e:
                logger.warning(
                    'Error occurred when parsing TC({},{}) due to {}'.format(
                        header['service_type'], header['service_subtype'],
                        str(e)))
        info = STIX_IDB.get_telecommand_info(header)
        if not info:
            logger.error(
                'Failed to retrieve telecommand information from the database')
            return stix_global.HEADER_KEY_ERROR, header

        header.update({
            'descr': info['CCF_DESCR'],
            'DESCR2': info['CCF_DESCR2'],
            'SPID': '',
            'name': info['CCF_CNAME'],
            'TMTC': 'TC',
            'SCET': 0
        })
        if status == stix_global.OK:
            try:
                header['ack_desc'] = stix_header.ACK_MAPPING[header['ack']]
            except KeyError:
                logger.error('Error occurs when  getting ack value')
                status = stix_global.HEADER_KEY_ERROR
        return status, header

    def parse_service_20(self, parameters):
        """
          The detailed structure of S20 is not defined in IDB
        """
        try:
            param_hb = parameters[9]
            #platform data parser to be added
        except IndexError:
            logger.warning(
                'Failed to extract STIX instrument data from a S20 packet')
            return
        param_obj = Parameter(param_hb)
        raw_bin = b'\x04' + param_obj.raw
        hb_param_children = self.parse_fixed_telemetry_packet(raw_bin, 54103)
        #parse them as for HK4
        #set list as its children
        param_obj.set_children(hb_param_children)

    def parse_binary(self, buf, i=0, header_auxiliary=None):
        """
        Inputs:
            buffer, i.e., the input binary array
            i: starting offset
            header_auxiliary: information attached to header, it must a dictionary
        Returns:
            decoded packets in python list if ret_packets=True
        """
        packets = []
        length = len(buf)
        self.inc_counter('total_length', length)
        if i >= length:
            return []
        raw_binary = b''
        while i < length and not self.stop_parsing:
            packet = None
            STIX_DECOMPRESSOR.reset()
            if buf[i] in stix_header.TM_HEADER_FIRST_BYTE:
                status, i, header_raw = get_from_bytearray(buf, i, 16)
                if status == stix_global.EOF:
                    break
                header_status, header = self.parse_telemetry_header(header_raw)
                if header_status != stix_global.OK:
                    logger.warning('Bad header at {}, code {} '.format(
                        i, header_status))
                    self.inc_counter('num_bad_headers')
                    continue

                self.inc_counter('num_tm')

                data_field_length = header['length'] - 9
                status, i, data_field_raw = get_from_bytearray(
                    buf, i, data_field_length)
                if status == stix_global.EOF:
                    logger.warning(
                        "Incomplete packet, the last {} bytes  were not parsed"
                        .format(len(data_field_raw)))
                    break
                ret = self.parse_data_field_header(header, data_field_raw,
                                                   data_field_length)
                if ret != stix_global.OK:
                    logger.warning(
                        'Missing information in the IDB to decoded the data starting at {} '
                        .format(i))
                    continue
                tpsd = header['TPSD']
                spid = header['SPID']
                self.inc_counter('spid', spid)
                if self.selected_services:
                    if header['service_type'] not in self.selected_services:
                        self.inc_counter('num_filtered')
                        continue
                if self.selected_spids:
                    if spid not in self.selected_spids:
                        self.inc_counter('num_filtered')
                        continue

                STIX_DECOMPRESSOR.init(spid)

                parameters = None
                if tpsd == -1:
                    #it is a fixed length telemetry packet
                    parameters = self.parse_fixed_telemetry_packet(
                        data_field_raw, spid)
                else:
                    # variable length telemetry packet
                    num_read, parameters, status = self.vp_tm_parser.parse(
                        data_field_raw, spid)
                    if num_read != data_field_length:
                        logger.warning(
                            ' Packet (SPID {}) data field size: {}B, actual read: {}B'
                            .format(spid, data_field_length, num_read))

                if header['service_type'] == 5 or (
                        header['service_type'] == 1
                        and header['service_subtype'] in [2, 8]):
                    self.stix_alerts.append(header)

                packet = {'header': header, 'parameters': parameters}
                self.inc_counter('num_tm_parsed')
                raw_binary = header_raw + data_field_raw
                #if self.store_binary:
                #    packet['bin'] = raw_binary

                #unique id

            elif buf[i] in stix_header.TC_HEADER_FIRST_BYTE:

                if len(buf) - i < 10:
                    logger.warning(
                        "Incomplete packet. The last {} bytes  were not parsed"
                        .format(len(buf)))
                    break
                header_status, header = self.parse_telecommand_header(buf, i)
                i += 10
                if header_status != stix_global.OK:
                    self.inc_counter('num_bad_headers')
                    logger.warning(
                        "Invalid telecommand header. ERROR code: {}, Current cursor at {} "
                        .format(header_status, i - 10))
                    continue
                self.inc_counter('num_tc')
                data_field_length = header['length'] + 1 - 4
                status, i, data_field_raw = get_from_bytearray(
                    buf, i, data_field_length)
                if status == stix_global.EOF:
                    break
                if self.selected_services:
                    if header['service_type'] not in self.selected_services:
                        self.inc_counter('num_filtered')
                        continue
                telecommand_name = header['name']
                num_read, parameters, status = self.tc_parser.parse(
                    telecommand_name, data_field_raw)
                if num_read != data_field_length - 2:  #the last two bytes is CRC
                    logger.warning(
                        ' TC {} data field size: {}B, actual read: {}B'.format(
                            telecommand_name, data_field_length, num_read))

                if telecommand_name == 'ZIX20128' and parameters:
                    #S20 detailed structure  not defined in ICD
                    if self.S20_excluded:
                        continue
                    self.parse_service_20(parameters)

                if header[
                        'service_type'] == 5 and header['service_subtype'] > 1:
                    self.instrument_message.append(header)
                    #used for instrument health tracking
                packet = {'header': header, 'parameters': parameters}
                self.inc_counter('num_tc_parsed')
                raw_binary = buf[i:i + 10] + data_field_raw
                #if self.store_binary:
                #    packet['bin'] = raw_binary

            else:
                old_i = i
                logger.warning('Unrecognized byte: {} at Pos {}'.format(
                    hex(buf[i]), i))
                i = find_next_header(buf, i)
                if i == stix_global.EOF:
                    break
                self.inc_counter('num_bad_bytes', i - old_i)
                logger.warning(
                    'New header found at Pos {}, {} bytes ignored!'.format(
                        i, i - old_i))

            if packet:
                packet['hash'] = hashlib.shake_256(raw_binary).hexdigest(8)
                if header_auxiliary:
                    for key, val in header_auxiliary.items():
                        packet['header'][key] = val
                else:
                    self.attach_timestamps(packet)

                if self.store_packet_enabled:
                    packets.append(packet)
                if self.packet_writer:
                    self.packet_writer.write_one(packet)

            logger.progress(i, length)

        return packets

    #def parse_live_hex_stream(self, raw):
    #    #used to parse live hex stream from TSC
    #    self.is_live_hex_stream = True
    #    #self.set_store_binary_enabled(False)
    #    self.set_packet_buffer_enabled(False)
    #    self.parse_hex(raw)

    def parse_hex(self, raw_hex):
        hex_string = re.sub(r"\s+", "", raw_hex)
        try:
            raw = binascii.unhexlify(hex_string)
        except Exception as e:
            logger.error(str(e))
            return []
        return self.parse_binary(raw)

    def parse_hex_file(self, filename):
        with open(filename, 'r') as f:
            raw_hex = f.read()
            return self.parse_hex(raw_hex)
        return None

    def parse_moc_ascii(self, filename):
        packets = []
        logger.set_progress_enabled(False)
        total_num_lines = sum(1 for line in open(filename))
        with open(filename) as filein:
            logger.info('Reading packets from the file {}'.format(filename))
            idx = 0
            for line in filein:
                try:
                    [self.receipt_utc, data_hex] = line.strip().split()
                    data_binary = binascii.unhexlify(data_hex)
                except Exception as e:
                    logger.error(str(e))
                    continue

                packet = self.parse_binary(data_binary)
                if packet:
                    packets.extend(packet)
                logger.set_progress_enabled(True)
                logger.progress(idx, total_num_lines)
                logger.set_progress_enabled(False)
                #not to show progress bar for each packet
                idx += 1
        return packets

    def attach_timestamps(self, packet):
        pkt_header = packet['header']

        pkt_header['UTC'] = ''
        pkt_header['unix_time'] = 0
        #if self.is_live_hex_stream:
        #    pkt_header['unix_time'] = stix_datetime.get_now('unix')
        #    pkt_header['UTC'] = stix_datetime.get_now('utc')
        #    return
        use_receipt_time = False
        if self.receipt_utc:
            #GU data
            try:
                dt = dtparser.parse(self.receipt_utc)
                pkt_header['UTC'] = dt
                pkt_header['unix_time'] = dt.timestamp()
                use_receipt_time = True
            except ValueError:
                logger.warning('Failed to parse timestamp: {}'.format(
                    self.receipt_utc))

        if pkt_header['TMTC'] == 'TM':
            coarse = pkt_header['coarse_time']
            fine = pkt_header['fine_time']
            pkt_header['obt_utc'] = stix_datetime.scet2utc(coarse, fine)
            if not use_receipt_time:
                pkt_header['unix_time'] = stix_datetime.scet2unix(coarse, fine)
                pkt_header['UTC'] = stix_datetime.unix2utc(
                    pkt_header['unix_time'])

    def parse_telecommand_report_xml(self, filename):
        #parse TC history
        #raw data should be included in the xml file
        packets = []
        results = []
        logger.set_progress_enabled(False)
        state_names = [
            'ReleaseState', 'GroundState', 'UplinkState', 'OnBoardState',
            'OnBoardAccState', 'OnBoardAccPBState', 'ExecCompPBState'
        ]
        with open(filename) as f:
            doc = xmltodict.parse(f.read())
            for item in doc['ns2:ResponsePart']['Response'][
                    'PktTcReportResponse']['PktTcReportList'][
                        'PktTcReportListElement']:
                if 'ZIX' not in item['CommandName'] or (
                        'RawBodyData' not in item):
                    #skip platform commands
                    continue
                state = ''.join([item[e][0] for e in state_names])
                header_auxiliary = {
                    'UTC': item['ExecutionTime'],
                    'unix_time': stix_datetime.utc2unix(item['ExecutionTime']),
                    'release_time': item['ReleaseTime'],
                    'uplink_time': item['UplinkTime'],
                    'execution_time': item['ExecutionTime'],
                    'sequence_name': item['SequenceName'],
                    'state': state,
                }
                raw = item['RawBodyData']
                packets.append({
                    'header_auxiliary': header_auxiliary,
                    'raw': raw
                })
        num = len(packets)
        for i, packet in enumerate(packets):
            data_hex = packet['raw']
            data = binascii.unhexlify(data_hex)
            result = self.parse_binary(data, 0, packet['header_auxiliary'])
            logger.set_progress_enabled(True)
            logger.progress(i, num)
            logger.set_progress_enabled(False)
            if not result:
                continue
            results.extend(result)
        return results

    def parse_telemetry_xml(self, raw_filename):
        packets = []
        results = []
        logger.set_progress_enabled(False)
        with open(raw_filename) as filein:
            logger.info('Parsing {}'.format(raw_filename))
            doc = xmltodict.parse(filein.read())
            for e in doc['ns2:ResponsePart']['Response']['PktRawResponse'][
                    'PktRawResponseElement']:
                packet = {'id': e['@packetID'], 'raw': e['Packet']}
                packets.append(packet)
        num = len(packets)
        for i, packet in enumerate(packets):
            data_hex = packet['raw']
            data_binary = binascii.unhexlify(data_hex)
            data = data_binary[76:]
            result = self.parse_binary(data)
            logger.set_progress_enabled(True)
            logger.progress(i, num)
            logger.set_progress_enabled(False)
            if not result:
                continue
            results.extend(result)
        return results

    def parse_file(self, raw_filename, file_type=None, clear=True):
        packets = []
        if clear:
            self.reset_counter()

        if self.packet_writer:
            self.packet_writer.set_filename(raw_filename)

        logger.info('Processing file: {}'.format(raw_filename))
        self.raw_filename = raw_filename
        self.in_filesize = os.path.getsize(raw_filename)

        if not file_type:
            file_type = detect_filetype(raw_filename)

        if file_type == 'bin':
            with open(raw_filename, 'rb') as in_file:
                data = in_file.read()
                packets = self.parse_binary(data)
        elif file_type == 'ascii':
            packets = self.parse_moc_ascii(raw_filename)
        elif file_type == 'xml.tm.raw':
            packets = self.parse_telemetry_xml(raw_filename)
        elif file_type == 'xml.tc.hist':
            packets = self.parse_telecommand_report_xml(raw_filename)
        elif file_type == 'hex':
            packets = self.parse_hex_file(raw_filename)
        else:
            logger.error('{} has unknown input file type'.format(raw_filename))
            return []

        summary = self.get_summary()
        logger.print_summary(summary)
        if self.packet_writer:
            self.packet_writer.set_summary(summary)
        return packets

    #def set_pickle_writer(self, out_filename, comment=''):
    #    self.packet_writer = stix_writer.StixPickleWriter(out_filename)
    #    idb_version = STIX_IDB.get_idb_version()
    #    self.packet_writer.register_run(self.raw_filename, self.raw_filename,
    #                                    comment, idb_version)

    def config_mongodb(self,
                           server,
                           port,
                           username,
                           password,
                           comment='',
                           raw_filename='',
                           instrument=''):
        #instrument: GU or PFM
        #server, port, username and password are required by MongoDB
        self.packet_writer = stix_writer.StixMongoDB(
            server, port, username, password)
        idb_version = STIX_IDB.get_idb_version()
        self.raw_filename = raw_filename
        self.packet_writer.register_run(raw_filename, self.in_filesize,
                                        comment, idb_version, instrument)

    def is_processed(self, filename):
        return self.packet_writer.is_processed(filename)

    def done(self):
        '''
            write a summary to MongoDB
            return the summary if success otherwise return None 

        '''
        if self.packet_writer:
            return self.packet_writer.close()
        return None

    def set_verbose_level(self, verbose_level):
        logger.set_level(verbose_level)

    def set_store_binary_enabled(self, stat):
        self.store_binary = stat

    def set_progress_bar_enabed(self, value):
        logger.set_progress_enabled(value)
