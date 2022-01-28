#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
# @title        : xml_to_bin.py
# @description  : read hex string from clipboard and write it to a binary file
# @date         : March. 28, 2019
import sys
import binascii
import tkinter as tk
import re
from io import BytesIO


def write_clipboard(fname):
    with open(fname, 'wb') as fout:
        root = tk.Tk()
        root.withdraw()
        raw_hex = root.clipboard_get()
        data_hex = re.sub(r"\s+", "", raw_hex)
        if len(data_hex) < 14:
            print('header:%s' % data_hex[:10])
            print('data invalid')
        else:
            data_binary = binascii.unhexlify(data_hex)
            fout.write(data_binary)
            print('data length:%d' % len(data_binary))
            print('written to:%s' % fname)


if __name__ == '__main__':
    fname = 'clipboard.dat'
    if len(sys.argv) == 2:
        fname = sys.argv[1]
    write_clipboard(fname)
