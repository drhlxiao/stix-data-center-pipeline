#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : parameter_description.py
# @description  : print the description  of a parameter
# @author       : Hualin Xiao
# @date         : March. 15, 2019
#

from __future__ import (absolute_import, unicode_literals)
import sys
from core import idb

_idb = idb._idb


def main():
    if len(sys.argv) == 1:
        print('Print the description of a parameter')
        print('Usage:')
        print('parameter_description  <parameter name>')
    else:
        parameter = sys.argv[1]
        print('\n\n')
        print(parameter + ':')
        desc = _idb.get_scos_description(parameter)
        if desc:
            print(str(desc))
        else:
            desc = _idb.get_PCF_description(parameter)
            if desc:
                print(str(desc[0]))


if __name__ == '__main__':
    main()
