#!/usr/bin/python3
import os
import sys
import binascii
import socket
import time
from datetime import datetime
from parser.core import writer
from PyQt5 import uic, QtWidgets, QtCore, QtGui
from PyQt5.QtCore import QThread, pyqtSignal, QTimer
from parser.core import parser
from parser.core import logger

class StixSocketPacketReceiver(QThread):
    """
    QThread to receive packets via socket
    """
    packetArrival = pyqtSignal(str)
    def __init__(self):
        super(StixSocketPacketReceiver, self).__init__()
        self.working = True
        self.port = 9000
        self.host = 'localost'

        #self.tctm_parser.set_report_progress_enabled(False)
        self.s = None

    def connect(self, host, port):
        self.working = True
        self.port = port
        self.host = host
        try:
            self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.s.connect((self.host, self.port))
        except Exception as e:
            print(str(e))

    def run(self):
        i=0
        if not self.s:
            return
        while True:
            buf = b''
            while True:
                data = self.s.recv(1)
                buf += data
                if data == b'>':
                    if buf.endswith(b'<-->'):
                        break
            data2 = buf.split()
            if buf[0:9] == 'TM_PACKET'.encode():
                data_hex = data2[-1][0:-4]
                hex_str=data_hex.decode()
                self.packetArrival.emit(hex_str)

class Main(QtCore.QObject):

    def __init__(self, parent=None):            
        super(Main, self).__init__(parent)   

        self.parser= parser.StixTCTMParser()
        self.parser.set_store_packet_enabled(False)
        self.parser.set_store_binary_enabled(False)

    def start_run(self,comment=''):
        self.thread= StixSocketPacketReceiver()
        self.thread.packetArrival.connect(self.processPacket)
        self.thread.connect('localhost',9001)
        self.thread.start()
        self.num_packets=0
        if comment:
            comment+=','
        comment+='live stream starts at: '+ str(datetime.now())
        self.parser.set_MongoDB_writer('localhost',27017, '','',comment)

    def processPacket(self, text):
        self.num_packets+=1
        print('Packet received:',self.num_packets)
        packets = self.parser.parse_live_hex_stream(text)
    def __del__(self):
        print("Total number of packet received:")
        print(self.num_packets)
        print('Done')
        self.parser.done()

if __name__ == "__main__":
    import sys
    import signal
    comment=''
    app = QtCore.QCoreApplication(sys.argv)
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    if len(sys.argv)>1:
        comment=sys.argv[1]
    a = Main()
    a.start_run(comment)
    app.exec()

