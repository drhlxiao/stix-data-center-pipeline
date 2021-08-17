#!/usr/bin/python3
# @author       : Hualin Xiao
# @date         : May. 11, 2021

import sys
import os
import time
import threading
from stix.pipeline import parser_pipeline as pd
from stix.pipeline import goes_downloader as gd

PARSER_SLEEP_TIME=60
GOES_TIME_LOOP=24*3600

def parser_loop():
    while True:
        num_new=pd.main()
        if num_new>0:
            gd.main()
        time.sleep(PARSER_SLEEP_TIME)

def goes_loop():
    while True:
        try:
            gd.main()
        except Exception as e:
            print(e)
        time.sleep(GOES_TIME_LOOP)



if __name__=='__main__':
    parser_thread= threading.Thread(target=parser_loop)
    parser_thread.start()
    goes_thread= threading.Thread(target=goes_loop)
    goes_thread.start()
