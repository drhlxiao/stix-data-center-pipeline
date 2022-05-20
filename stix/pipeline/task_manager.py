#!/usr/bin/env python 
#test.py

import multiprocessing
import time
from stix.core import logger

logger = logger.get_logger()

def run_on_background(func, name="process", args=(), wait=3600*3):
    logger.info(f"Running {name} on background!")
    p = multiprocessing.Process(target=func, name=name, args=args)
    p.start()
    p.join(wait)
    if p.is_alive():
        #kill it after time out
        logger.info(f"Timeout Terminating {name} !")
        p.terminate()
