#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : IOR2tex.py
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
import asw_param
MAX_PARAM_DISPLAY = 10
import re


def get_subsystem_name(system, subID):
    if system == 1 or system == "IDPU":
        return "Watchdog"
    elif system == 2 or system == "DETECTORS":
        mmap = {
            1: 'detector quadrant 1',
            2: 'detector quadrant 2',
            4: 'detector quadrant 3',
            8: 'detector quadrant 4'
        }
        return mmap[subID]
    elif system == 3 or system == "ASPECT_SYSTEM":
        if subID == 1:
            return "Aspect system A"
        return "Aspect system B"
    elif system == 4 or system == "ATTENUATOR":
        if subID == 1:
            return "Motor 1"
        if subID == 2:
            return "Motor 2"
        if subID == 4:
            return "Override Motor 1"
        if subID == 8:
            return "Override Motor 2"
    elif system == 5 or system == "PSU":
        if (subID == 4):
            return "LV power supply"
        if (subID == 1):
            return "HV power supply 1"
        if (subID == 2):
            return "HV power supply 2"

    return ""


def tex_escape(text):
    """
        :param text: a plain text message
        :return: the message escaped to appear correctly in LaTeX
    """
    conv = {
        '&': r'\&',
        #'%': r'\%',
        #'$': r'\$',
        #'#': r'\#',
        '_': r'\_',
        #'{': r'\{',
        #'}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\^{}',
        '<': r'\textless{}',
        '>': r'\textgreater{}',
    }
    regex = re.compile('|'.join(
        re.escape(str(key))
        for key in sorted(conv.keys(), key=lambda item: -len(item))))
    return regex.sub(lambda match: conv[match.group()], text)


#    if param_name=='PIX00016':
#        return asw_param.get_name(value)
#    elif param_name=='PIX00083':
#        return ' ({})'.format(get_subsystem_name(


def get_parameters(parameters):
    if not parameters:
        return ''
    content = '\\begin{itemize}\n'
    i = 0
    for param in parameters:
        eng = param[2]
        content += '\\item  {} ({}) = {} ({})\n'.format(
            param[4], param[2], param[1], eng)
        i += 1
        if i == 10:
            break

    content += '\\end{itemize}'
    return content


def read_PDOR(filename):
    reader = ior2dict.IORReader()
    result = reader.parse(filename)
    occurrences = result['occurrences']
    content = 'The PDOR \\pdor({}) contains steps as follow:\n'.format(
        filename)
    content += '\\begin{itemize}\n'
    for occ in occurrences:
        desc = occ['desc'][0]
        name = occ['name']
        if not desc:
            desc = occ['desc'][1]

        parames = get_parameters(occ['parameters'])

        item = ' '.join(
            ['\\item Execute \\textbf{', name, '}', desc, '\n', parames, '\n'])
        content += tex_escape(item)

    content += '\\end{itemize}\n'
    return content


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('convert PDOR to tex script')
        print('./ior2doc pdor')
    else:
        tex_filename = sys.argv[1].replace('SOL', 'tex')
        with open(tex_filename, 'w') as f:
            content = read_PDOR(sys.argv[1])
            f.write(content)
            print('Output:', tex_filename)
