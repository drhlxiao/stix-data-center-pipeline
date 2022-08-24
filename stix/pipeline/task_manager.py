#!/usr/bin/env python 
#test.py

import multiprocessing
import time

from stix.core import logger
from threading import Timer
logger = logger.get_logger()
def kill_me(p, name):
    logger.info(f"Timeout! Terminating {name}  !")
    try:
        if p.is_alive():
            logger.info(f"Killing it!")
            p.terminate()
    except Exception:
        pass

def run_on_background(func, name="process", args=(), wait=3600*3):
    logger.info(f"Running {name} on background!")
    p = multiprocessing.Process(target=func, name=name, args=args)
    p.start()
    t=Timer(wait, kill_me, [p, name])
    t.start()

