#!/usr/bin/env python 
#test.py

import multiprocessing
import time
from stix.core import logger

logger = logger.get_logger()

def run_on_background(func, name='processing', args, wait=3600*3):
    logger.info(f"Running {name} on background!")
    
    p = multiprocessing.Process(target=func, name=name, args)
    p.start()
    time.sleep(wait)
    logger.info(f"Terminating {name} !")
    p.terminate()
    p.join()
