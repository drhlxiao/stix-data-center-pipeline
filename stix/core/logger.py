#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# @title        : logger.py
# @description  : logger

import sys
from datetime import datetime

DEBUG = 1
PROGRESS = 2
INFO = 3
WARNING = 4
ERROR = 5
CRITICAL = 6


class StixLogger(object):
    __instance = None

    @staticmethod
    def get_instance(filename=None, level=4):
        if not StixLogger.__instance:
            StixLogger(filename, level)
        return StixLogger.__instance

    #singleton

    def __init__(self, filename=None, level=10):

        if StixLogger.__instance:
            raise Exception('Logger already initialized')
        else:
            StixLogger.__instance = self
        self.logfile = None
        self.filename = filename
        self.set_logger(filename, level)
        self.debug = True

        self.progress_enabled = True

        self.last_progress = 0

        self.signal_handler = None
        self.progress_bar_last_num = 0

    def set_progress_enabled(self, status):
        self.progress_enabled = status

    def set_signal_handlers(self, handler):
        self.signal_handler = handler

    def get_now(self):
        now = datetime.now()
        return now.strftime("%Y-%m-%d %H:%M:%S")

    def emit(self, msg, level):
        if not self.signal_handler:
            return
        self.signal_handler[level].emit(msg)

    def update_console_progress(self, progress):
        num = int(progress / 2.)
        if num > self.progress_bar_last_num:
            bar = '#' * num + ' ' * (50 - num)
            #print is slow
            print('\r[{0}] {1}%'.format(bar, progress), end="", flush=True)
        self.progress_bar_last_num = num

    def get_log_filename(self):
        return self.filename

    def enable_debugging(self):
        self.debug = True

    def set_logger(self, filename=None, level=3):
        if self.logfile:
            self.logfile.close()
            self.logfile = None
        self.filename = filename
        self.level = level
        if filename:
            try:
                self.logfile = open(filename, 'w+')
            except IOError:
                print('Can not open log file {}'.format(filename))

    def set_level(self, level):
        self.level = level

    def write(self, msg, level=INFO):
        if self.debug:
            print(msg)
        elif self.signal_handler:
            self.emit(msg, level)
        elif self.logfile:
            #if level == PROGRESS:
            #    msg = '{}% processed.'.format(msg)
            self.logfile.write('{}\n'.format(msg))
        else:
            print(msg)

    def critical(self, msg):
        self.write(('[INFO {}] : {}'.format(self.get_now(), msg)), CRITICAL)

    def error(self, msg):
        self.write(('[ERROR {}] : {}'.format(self.get_now(), msg)), ERROR)
    def warn(self,msg):
        self.warning(msg)
    def alert(self,msg):
        self.warning(msg)

    def warning(self, msg):
        if self.level < WARNING:
            return
        self.write(('[WARN {}] : {}'.format(self.get_now(), msg)), WARNING)

    def info(self, msg):
        if self.level < INFO:
            return
        self.write(('[INFO {}] : {}'.format(self.get_now(), msg)), INFO)

    def debug(self, msg):
        if self.level < DEBUG:
            return
        self.write(msg)

    def progress(self, i, total):
        if self.progress_enabled:
            current = int(100. * i / total)
            if current > self.last_progress:
                if self.logfile:
                    self.write(current, INFO)
                elif self.signal_handler:
                    self.emit(current, PROGRESS)
                else:
                    self.update_console_progress(current)

            self.last_progress = current

    def print_summary(self, summary):
        self.critical(
            'Size: {} bytes (bad:{});'
            ' Nb. of packets: {} ('
            'TM: {}, TC:{}, Filtered: {}); Parsed {} (TM:{},TC:{}); Bad headers:{} .'
            .format(
                summary['total_length'], summary['num_bad_bytes'],
                summary['num_tm'] + summary['num_tc'] +
                summary['num_filtered'], summary['num_tm'], summary['num_tc'],
                summary['num_filtered'],
                summary['num_tm_parsed'] + summary['num_tc_parsed'],
                summary['num_tm_parsed'], summary['num_tc_parsed'],
                summary['num_bad_headers']))


def get_logger(filename=None):
    return StixLogger.get_instance(filename)
