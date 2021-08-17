#!/usr/bin/python3
import time
import os
import sys
import re
import ior2dict

from datetime import datetime

DESCRIPTION = """
A tool to validate STIX PDOR Printouts.
Usage:
    printout_validate PDOR_DIRECTORY  [PRINTOUT_DIRECTORY]
"""

import sys
import os
from pprint import pprint
import csv

html_template = '''
<!doctype html>
<html lang="en">
<head>
  <title>Printout  validation</title>
  <style>
  .warning
  {{
    color : #FFC300;
  }}
  .error
  {{
    color : red;
  }}
  .ok
  {{
    color : green;
  }}
  .row-0
  {{
    font-weight:bold;
   }}
   .col-1 
   {{
    font-weight:bold;
   }}
</style>

</head>
<body>
{}
<h3 class="{}">{}</h3>
<p>Created at {} by printout validation tool</p>

</body>
</html>
'''
num_errors = 0
num_warnings = 0

ior_reader = ior2dict.IORReader()


#pprint(CSS)
def error(msg):
    global num_errors
    num_errors += 1
    print('[ERROR]:{}'.format(msg))


def info(msg):
    print('[INFO]:{}'.format(msg))


def warn(msg):
    global num_warnings
    num_warnings += 1
    print('[WARN]:{}'.format(msg))


def success(msg):
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


def are_parameters_equal(param1, param2):
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


def check_parameters(params1, params2, rows):
    mismatches = []

    for p1 in params1:
        if not p1:
            continue
        row = [''] * 8
        row[2] = p1[4]
        row[3] = p1[1]

        match = False
        for p2 in params2:
            row[4] = p2[0]
            row[5] = p2[3]
            if are_parameters_equal(p1, p2):
                match = True
                break
        if not match:
            row[7] += '<span class="error">invalid</span>'
            mismatches.append(p1[4])
        else:
            row[7] += '<span class="ok">valid</span>'
        rows.append(row)

    return mismatches


def to_seconds(action_time):
    if action_time in ['ASAP', 'asap']:
        return 0

    is_absolute = re.findall(r'\d{4}-\d{2}-\d{2}', action_time)
    if is_absolute:
        return dtparser.parse(action_time).timestamp()
    else:
        if re.findall(r'\d{2}:\d{2}:\d{2}', action_time):
            time_fields = action_time.split(':')
            return int(time_fields[0]) * 3600 + int(time_fields[1]) * 60 * int(
                time_fields[2])
        if re.findall(r'\d{2}.\d{2}.\d{2}', action_time):
            time_fields = action_time.split('.')
            return int(time_fields[0]) * 3600 + int(time_fields[1]) * 60 * int(
                time_fields[2])
    return -1


def check_timestamps(occ_time, pro_time):
    if occ_time == pro_time:
        return True
    occ_time_seconds = to_seconds(occ_time)
    pro_time_seconds = to_seconds(pro_time)

    if occ_time_seconds >= 0 and pro_time_seconds >= 0 and occ_time_seconds == pro_time_seconds:
        return True
    return False


def check_occurrence(i, occ, printouts, row):
    pro = None
    info('checking #{}: {}'.format(i, occ['name']))
    while True:
        pro = printouts.pop(0)
        if 'ZIX' in pro['name'] or 'AIX' in pro['seq']:
            #skip platform telecommands
            break
        else:
            info('Ignored TC in printouts:{}/{}'.format(
                pro['name'], pro['seq']))

    row[4] = pro['name']
    info('checking #{}: {}'.format(i, pro['name']))

    status = 'ok'
    if not check_timestamps(occ['actionTime'], pro['duration']):
        status = 'warning'
        warn('TC {} action time {} <--> {} are not the same'.format(
            occ['name'], occ['actionTime'], pro['duration']))
    row[6] = '<span class="{};" >{}</span>'.format(status, pro['duration'])
    param_rows = []

    if pro['name'] == occ['name'] or pro['seq'] == occ['name']:
        row[7] = '<span class="ok" >valid</span>'
        paramter_mismatches = check_parameters(occ['parameters'],
                                               pro['parameters'], param_rows)
        if paramter_mismatches:
            error('Mismatched parameters in {}: {}'.format(
                occ['name'], str(paramter_mismatches)))
        else:
            success('{} validated'.format(occ['name']))
    else:
        error(
            'In PDOR cmd {} is not the same as in printout:  {} - {}  '.format(
                occ['name'], pro['name'], pro['seq']))
        row[7] = '<span class="error" >Two telecommands are not the same</span>'
    return param_rows


def validate_pdor_printout(pdor_filename, printout_filename):
    pdors = ior_reader.parse(pdor_filename, no_seq=True)
    printouts = read_printouts(printout_filename)

    content = '<br><br><h3> Diff. {} and {}</h3>\n'.format(
        printout_filename, pdor_filename)
    table_template = '<table border="1">{}</table>\n'
    table = '<tr class="row-0"><td>#</td><td colspan="3">{}</td> <td colspan="4">{}</td></tr>'.format(
        printout_filename, pdor_filename)

    table_rows = [[
        '', 'TC name', 'description', 'action time/value', 'TC name',
        'description', 'action_time/value', 'checks'
    ]]
    param_rows = []
    for i, occ in enumerate(pdors['occurrences']):
        row = [i, occ['name'], occ['desc'], occ['actionTime'], '', '', '', '']
        if 'ZIX' not in occ['name'] and 'AIX' not in occ['name']:
            success('Platform command {} ignored'.format(occ['name']))
            row[7] = 'Platform command ignored'
        else:
            param_rows = check_occurrence(i, occ, printouts, row)
        table_rows.append(row)
        if param_rows:
            table_rows.extend(param_rows)

    for irow, tr in enumerate(table_rows):
        table += '<tr class="row-{}">'.format(irow)
        for icol, col in enumerate(tr):
            table += '<td class="col-{}">{}</td>'.format(icol, col)
        table += '</tr>'

    return content + table_template.format(table)


def validate(pdor_path, printout_path=None):
    cfiles = []
    if not printout_path:
        printout_path = pdor_path
    html_filename = os.path.join(printout_path, 'printouts_check.html')
    body = ''
    report = ''
    with open(html_filename, 'w') as html:
        for root, dirs, files in os.walk(printout_path):
            for file in files:
                if file.endswith('_F.csv') and file.startswith('PDOR'):
                    printout_filename = os.path.join(root, file)
                    pdor_filename = os.path.join(
                        pdor_path,
                        os.path.basename(printout_filename[:-6] + '.SOL'))
                    info('check file:{} <-> {}'.format(printout_filename,
                                                       pdor_filename))
                    if os.path.exists(pdor_filename):
                        report += validate_pdor_printout(
                            pdor_filename, printout_filename)
                    else:
                        msg = 'PDOR: {} does not exist'.format(pdor_filename)
                        warn(msg)
                        report += '<h3 class="warning">{}</h3>'.format(msg)

        summary = 'Total number of warnings :{}<br>\n Total number of errors {}'.format(
            num_warnings, num_errors)
        info(summary)
        status = 'ok'
        if num_errors > 0:
            status = 'error'
        timestamp = datetime.now().isoformat()
        body = html_template.format(report, status, summary, timestamp)
        html.write(body)
    info('check report:{}'.format(html_filename))


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Usage: ./validate_printouts PDOR_PATH  [PRINTOUT_PATH]')
    elif len(sys.argv) == 2:
        validate(sys.argv[1])
    else:
        print('#here')
        validate(sys.argv[1], sys.argv[2])
