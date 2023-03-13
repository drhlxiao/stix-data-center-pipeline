#!/usr/bin/python3
# @author       : Hualin Xiao
# @date         : Mar. 06, 2023
# @description:  manage a task


import time
import multiprocessing
from stix.core import logger
logger = logger.get_logger()



class TaskManager(object):
    def __init__(self, func,  sleep_time=120, timeout=7200, task_name=''):
        self.task_name=task_name
        self.sleep_time = sleep_time
        self.timeout = timeout
        self.func= func
    def start(self): 
        #start the ps
        logger.info(f'starting task {self.task_name}...')
        ps = multiprocessing.Process(target=self.func, name=self.task_name)
        ps.start()
        # Check the ps every 5 minutes
        ps.join(self.timeout)
        if ps.is_alive():
            logger.info("Time out..")
            ps.terminate()

        time.sleep(self.sleep_time)
        self.start()
        



