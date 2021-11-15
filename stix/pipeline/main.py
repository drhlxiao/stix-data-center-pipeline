#!/usr/bin/python3
# @author       : Hualin Xiao
# @date         : May. 11, 2021

import sys
import os
import time
import threading
from stix.pipeline import parser as pd
from stix.analysis import goes_downloader as gd

goes=gd.GOES()



def parser_loop(wait=60):
    while True:
        num_new=pd.main()
        time.sleep(wait)

def main_threading():
    parser_thread= threading.Thread(target=parser_loop)
    parser_thread.start()
    goes_thread= threading.Thread(target=goes.loop)
    goes_thread.start()

if __name__ == '__main__':
    main_threading()

