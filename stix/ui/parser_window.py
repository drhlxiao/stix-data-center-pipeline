#!/usr/bin/python3
# -*- encoding: utf-8 -*-

import os
import sys
import pickle
import time
import gzip
import binascii
import struct
import pprint
import socket
import signal
import webbrowser
from functools import partial
from stix.spice import time_utils as tu
from datetime import datetime
import numpy as np

from PyQt5 import  QtWidgets, QtCore, QtGui
from PyQt5.QtCore import QThread, pyqtSignal, QTimer
from PyQt5.QtChart import QChart, QChartView, QLineSeries, QValueAxis
from PyQt5.QtChart  import QBarSeries, QBarSet, QScatterSeries

from stix.core import parser
from stix.core import io as stixio
from stix.core import idb
from stix.core import mongo_db as mgdb
from stix.core import logger
from stix.ui import mainwindow
from stix.ui import mongo_dialog
from stix.ui import tsc_connection
from stix.ui import packet_filter
from stix.ui import plugin
from stix.ui import timestamp_convertor
from stix.ui import raw_viewer
#from stix.ui import console
from stix.core.datatypes import Parameter
from stix.core.datatypes import Packet

SELECTED_SERVICES = [1, 3, 5, 6, 9, 17, 20, 21, 22, 236, 237, 238, 239]

STIX_IDB = idb.idb()
MAX_NUM_PACKET_IN_BUFFER = 6000

LOGGER = logger.get_logger()


class ParserQThread(QThread):
    error = pyqtSignal(str)
    info = pyqtSignal(str)
    critical = pyqtSignal(str)
    warning = pyqtSignal(str)
    dataLoaded = pyqtSignal(list)
    progress = pyqtSignal(float)
    packetArrival = pyqtSignal(list)

    def __init__(self):
        super(ParserQThread, self).__init__()
        self.tctm_parser = parser.StixTCTMParser()
        self.tctm_parser.set_packet_buffer_enabled(True)
        self.tctm_parser.set_store_binary_enabled(True)

        handlers = {
            logger.INFO: self.info,
            logger.CRITICAL: self.critical,
            logger.WARNING: self.warning,
            logger.ERROR: self.error,
            logger.PROGRESS: self.progress
        }
        LOGGER.set_signal_handlers(handlers)

    def connectSignalSlots(self, slots):
        self.error.connect(slots['error'])
        self.info.connect(slots['info'])
        self.warning.connect(slots['warning'])
        self.dataLoaded.connect(slots['dataLoaded'])
        self.packetArrival.connect(slots['packetArrival'])
        self.critical.connect(slots['critical'])
        self.progress.connect(slots['progress'])


class StixSocketPacketReceiver(ParserQThread):
    """
    QThread to receive packets via socket
    """

    def __init__(self):
        super(StixSocketPacketReceiver, self).__init__()
        self.working = True
        self.port = 9000
        self.host = 'localost'

        LOGGER.set_progress_enabled(False)
        self.s = None

    def connect(self, host, port):
        self.working = True
        self.port = port
        self.host = host
        try:
            self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.s.connect((self.host, self.port))
            self.info.emit('Listening server {}:{}'.format(
                self.host, self.port))
        except Exception as e:
            self.error.emit(str(e))

    def run(self):
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
                data_binary = binascii.unhexlify(data_hex)
                packets = self.tctm_parser.parse_binary(data_binary)
                if packets:
                    packets[0]['header']['arrival'] = str(datetime.now())
                    self.packetArrival.emit(packets)


class StixSocketPacketServer(ParserQThread):
    """
    QThread to receive packets via socket
    """

    def __init__(self):
        super(StixSocketPacketServer, self).__init__()
        self.working = True
        self.sock = None
        self.current_id = 0
        self.packets = []
        self.running = False
        self.connection = None

    def setData(self, current_id, packets):
        self.current_id = current_id
        self.packets = packets

    def connect(self, host, port):
        if self.sock:
            return
        try:
            self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.sock.bind((host, port))
            self.sock.listen(1)
            self.info.emit('Packet server started: localhost:{}'.format(port))
        except Exception as e:
            self.error.emit(str(e))

    def run(self):
        EOF = b'EOF'
        if not self.sock:
            return
        if self.running:
            return
        try:
            while True:
                self.running = True
                self.connection, addr = self.sock.accept()
                self.info.emit('Connection from {}'.format(addr))
                control = self.connection.recv(64).decode('utf-8')

                if 'len' in control:
                    packet_length = len(self.packets)
                    data = pickle.dumps(packet_length) + EOF
                    length = len(data)
                    data_send = struct.pack('>I', length) + data
                    self.connection.sendall(data_send)
                else:
                    index = -1
                    slice_notation = ''
                    try:
                        index = int(control)
                    except ValueError:
                        if ':' in control:
                            slice_notation = control
                    pkts = []
                    try:
                        if index > -1:
                            pkts = [self.packets[index]]
                        if slice_notation:
                            pkts = eval('self.packets[{}]'.format(control))
                    except:
                        LOGGER.warn('Request from socket invalid')
                    data = pickle.dumps(pkts) + EOF
                    length = len(data)
                    size = struct.pack('>I', length)
                    data_send = size + data
                    self.info.emit('Sending {} kB to {} ...'.format(
                        round(length / 1024, 2), addr))
                    self.connection.sendall(data_send)
                    self.info.emit('{} kB sent.'.format(length / 1024))

                time.sleep(0.2)
        except Exception as e:
            self.error.emit(str(e))
            if self.connection:
                self.connection.close()
        finally:
            if self.connection:
                self.connection.close()


class StixFileReader(ParserQThread):
    """
    thread
    """

    def __init__(self):
        super(StixFileReader, self).__init__()
        self.data = []
        self.filename = None
        LOGGER.set_progress_enabled(True)

    def setFilename(self, filename):
        self.filename = filename
        self.data.clear()

    def setPacketFilter(self, selected_services, selected_SPID):
        self.tctm_parser.set_packet_filter(selected_services,
                                                selected_SPID)

    def stopParsing(self):
        self.tctm_parser.kill()

    def run(self):
        self.data = []
        filename = self.filename
        if filename.endswith('.pklz'):
            self.info.emit('Loading {} ...'.format(filename))
            f = gzip.open(filename, 'rb')
            self.data = pickle.load(f)['packet']
            f.close()
        elif filename.endswith('.pkl'):
            f = open(filename, 'rb')
            self.info.emit('Loading ...')
            self.data = pickle.load(f)['packet']
            f.close()
        else:
            self.data = self.tctm_parser.parse_file(filename)
        if self.data:
            self.dataLoaded.emit(self.data)


class StixHexStringParser(ParserQThread):
    def __init__(self):
        super(StixHexStringParser, self).__init__()
        self.data = []
        self.hex_string = []

    def setHex(self, hex_str):
        self.hex_string = hex_str
        self.data.clear()

    def run(self):
        if self.hex_string:
            self.data = self.tctm_parser.parse_hex(self.hex_string)
            if self.data:
                self.packetArrival.emit(self.data)


class Ui(mainwindow.Ui_MainWindow):
    def __init__(self, MainWindow):
        super(Ui, self).setupUi(MainWindow)
        self.MainWindow = MainWindow
        self.socketPacketReceiver = None

        self.timmer_is_on = False
        self.hexParser = StixHexStringParser()
        self.socketPacketReceiver = StixSocketPacketReceiver()
        self.socketPacketServer = StixSocketPacketServer()
        self.dataReader = StixFileReader()

        slots = {
            'info': self.onDataReaderInfo,
            'warning': self.onDataReaderWarning,
            'error': self.onDataReaderError,
            'critical': self.onDataReaderCritical,
            'dataLoaded': self.onDataReady,
            'packetArrival': self.onPacketArrival,
            'progress': self.onProgressUpdated
        }

        self.socketPacketReceiver.connectSignalSlots(slots)
        self.dataReader.connectSignalSlots(slots)
        self.hexParser.connectSignalSlots(slots)
        self.socketPacketServer.connectSignalSlots(slots)


        self.tabWidget.setCurrentIndex(0)
        self.actionExit.triggered.connect(self.close)
        self.actionPlot.setEnabled(False)

        self.actionNext.setIcon(self.style().standardIcon(
            QtWidgets.QStyle.SP_ArrowForward))
        self.actionPrevious.setIcon(self.style().standardIcon(
            QtWidgets.QStyle.SP_ArrowBack))
        self.actionOpen.setIcon(self.style().standardIcon(
            QtWidgets.QStyle.SP_DialogOpenButton))
        self.actionSave.setIcon(self.style().standardIcon(
            QtWidgets.QStyle.SP_DriveFDIcon))

        self.actionSave.triggered.connect(self.save)

        self.actionOpen.triggered.connect(self.getOpenFilename)
        self.filterPattern = None

        self.actionNext.triggered.connect(self.nextPacket)
        self.actionPrevious.triggered.connect(self.previousPacket)
        self.actionAbout.triggered.connect(self.about)

        self.actionPrevious.setEnabled(False)
        self.actionNext.setEnabled(False)
        self.actionSave.setEnabled(False)
        self.actionPlot.setEnabled(False)
        self.actionCopy.triggered.connect(self.onCopyTriggered)

        self.packetTreeWidget.currentItemChanged.connect(self.onPacketSelected)

        self.actionCopy.setEnabled(False)
        self.actionPaste.triggered.connect(self.onPasteTriggered)
        self.actionLog.triggered.connect(self.dockWidget.show)
        self.actionSetIDB.triggered.connect(self.onSetIDBClicked)
        self.plotButton.clicked.connect(
            partial(self.onPlotButtonClicked, None))

        #self.progressBar = QtWidgets.QProgressBar()
        #self.statusbar.addPermanentWidget(self.progressBar)
        self.progressDiag = None
        #self.progressBar.hide()

        self.actionPacketServer.triggered.connect(self.startPacketServer)

        self.exportButton.clicked.connect(self.onExportButtonClicked)
        self.actionPlot.triggered.connect(self.onPlotActionClicked)
        self.actionLoadMongodb.triggered.connect(self.onLoadMongoDBTriggered)
        self.actionConnectTSC.triggered.connect(self.onConnectTSCTriggered)
        self.actionPacketFilter.triggered.connect(self.filter)
        self.actionPlugins.triggered.connect(self.onPluginTriggered)
        self.actionOnlineHelp.triggered.connect(self.onOnlineHelpTriggered)
        self.actionViewBinary.triggered.connect(self.onViewBinaryTriggered)
        self.actionTimestampConvertor.triggered.connect(self.onTimestampConvertorTriggered)
        #self.actionPythonConsole.triggered.connect(self.startPythonConsole)
        self.autoUpdateButton.clicked.connect(
            self.onPlotAutoUpdateButtonClicked)

        self.packetTreeWidget.customContextMenuRequested.connect(
            self.packetTreeContextMenuEvent)

        #self.statusListWidget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #self.statusListWidget.customContextMenuRequested.connect(self.statusListContextMenuEvent)

        self.mdb = None

        self.current_row = 0
        self.data = []
        self.x = []
        self.y = []
        self.xlabel = 'x'
        self.ylabel = 'y'

        self.buttons_enabled = False

        self.chart = QChart()
        self.chart.layout().setContentsMargins(0, 0, 0, 0)
        self.chart.setBackgroundRoundness(0)
        self.savePlotButton.clicked.connect(self.savePlot)

        self.chartView = QChartView(self.chart)
        self.gridLayout.addWidget(self.chartView, 1, 0, 1, 15)
        self.selected_services = SELECTED_SERVICES
        self.selected_SPID = []
        self.selected_tmtc = 3

        # IDB location

        self.settings = QtCore.QSettings('FHNW', 'parser')
        self.idb_filename = self.settings.value('idb_filename', [], str)
        if self.idb_filename:
            STIX_IDB.reload(self.idb_filename)
        if not STIX_IDB.is_connected():
            self.showMessage('IDB has not been set!')
        else:
            idb_filename = STIX_IDB.get_idb_filename()
            self.showMessage('IDB loaded from : {} '.format(idb_filename), 1)
            if idb_filename != self.idb_filename:
                self.settings.setValue('idb_filename', idb_filename)
                self.idb_filename = idb_filename

    #def startPythonConsole(self):
    #    console.start({'packets': self.data})
    def close(self):
        self.MainWindow.close()
    def style(self):
        return self.MainWindow.style()

    def startPacketServer(self):
        host = 'localhost'
        port = 9096

        self.socketPacketServer.connect(host, port)
        self.socketPacketServer.setData(self.current_row, self.data)
        self.socketPacketServer.start()

        abspath = os.path.dirname(os.path.abspath(__file__))

        template = (
            "packets=req.request(query_str='len', host='{}',port={}, verbose_level=1)\n"
            "#a query_string can be \n"
            "#  -  a python slice notation, for example, ':' '0:-1', 3:-1\n"
            "#  -  'len',  to get the total number of packets,\n"
            "#  -   index ,  to get a packet of the given index"
            "#set verbose_level to 0, to suppress print output  ").format(
                abspath, host, port)
        cb = QtWidgets.QApplication.clipboard()
        cb.clear(mode=cb.Clipboard)
        cb.setText(template, mode=cb.Clipboard)
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText(
            "Packet server started and a template to request packet has been copied to your clipboard!"
        )
        retval = msg.exec_()

    def onPlotAutoUpdateButtonClicked(self):
        if not self.timmer_is_on:
            if not self.data:
                return
            num_packets = len(self.data)
            if num_packets > 200:
                packets = self.data[-200:-1]
            else:
                packets = self.data
            self.timer = QTimer()
            self.timer.timeout.connect(
                partial(self.onPlotButtonClicked, packets))
            self.timer.start(2000)
            self.timmer_is_on = True
            self.autoUpdateButton.setText('Stop Auto Update')
        else:
            if self.timer:
                self.timer.stop()
            self.timmer_is_on = False
            self.autoUpdateButton.setText('Start Auto Update')

    #def statusListContextMenuEvent(self,pos):
    #    menu = QtWidgets.QMenu()
    #    clearLogAction= menu.addAction('Empty log')
    #    clearLogAction.triggered.connect(self.clearLog)
    #def clearLog(self):
    #    self.statusListWidget.clear()

    def packetTreeContextMenuEvent(self, pos):
        menu = QtWidgets.QMenu()
        filterAction = menu.addAction('Filter')
        menu.addSeparator()
        rawDataAction = menu.addAction('Raw binary data')
        menu.addSeparator()
        copyPacketAction = menu.addAction('Copy packet')
        menu.addSeparator()
        deleteAllAction = menu.addAction('Delete all packets')
        self.current_row = self.packetTreeWidget.currentIndex().row()

        rawDataAction.triggered.connect(self.onViewBinaryTriggered)
        filterAction.triggered.connect(self.filter)
        copyPacketAction.triggered.connect(self.onCopyTriggered)
        deleteAllAction.triggered.connect(self.onDeleteAllTriggered)
        action = menu.exec_(self.packetTreeWidget.viewport().mapToGlobal(pos))

    def filter(self):
        text, okPressed = QtWidgets.QInputDialog.getText(
            None, "Packet filtering",
            "Filtering by SPID or description (! to exclude):",
            QtWidgets.QLineEdit.Normal, "")
        if okPressed:
            self.filterPattern = text
            self.addPacketsToView(self.data, True, show_stat=False)

    def onDeleteAllTriggered(self):
        self.data.clear()
        self.current_row = 0
        self.packetTreeWidget.clear()
        self.paramTreeWidget.clear()

    #def onPacketTreeItemDoubleClicked(self):
    #    self.onViewBinaryTriggered()

    def onViewBinaryTriggered(self):
        diag = QtWidgets.QDialog()
        diag_ui = raw_viewer.Ui_Dialog()
        diag_ui.setupUi(diag)
        if self.data:
            try:
                raw = self.data[self.current_row]['bin']
                header = self.data[self.current_row]['header']
                diag_ui.setPacketInfo('{}({},{})  {}'.format(
                    header['TMTC'], header['service_type'],
                    header['service_subtype'], header['descr']))
                diag_ui.displayRaw(raw)
            except (IndexError, KeyError):
                diag_ui.setText('Raw data not available.')
        diag.exec_()

    def onOnlineHelpTriggered(self):
        webbrowser.open(
            'https://github.com/i4Ds/STIX-python-data-parser', new=2)

    def onTimestampConvertorTriggered(self):
        diag = QtWidgets.QDialog()
        diag_ui = timestamp_convertor.Ui_Dialog()
        diag_ui.setupUi(diag)
        diag.exec_()


    def onPluginTriggered(self):
        self.plugin_location = self.settings.value('plugin_location', [], str)
        diag = QtWidgets.QDialog()
        diag_ui = plugin.Ui_Dialog()
        diag_ui.setupUi(diag)
        if self.plugin_location:
            diag_ui.setPluginLocation(self.plugin_location)

        diag_ui.setData(self.data, self.current_row)
        diag.exec_()
        location = diag_ui.getPluginLocation()
        if location != self.plugin_location:
            self.settings.setValue('plugin_location', location)

    def onPacketFilterTriggered(self):
        diag = QtWidgets.QDialog()
        diag_ui = packet_filter.Ui_Dialog()
        diag_ui.setupUi(diag)
        self.filterPattern = ''  #empty search string
        diag_ui.setSelectedServices(self.selected_services)

        diag_ui.buttonBox.accepted.connect(
            partial(self.applyServiceFilter, diag_ui))
        diag.exec_()

    def applyServiceFilter(self, diag_ui):
        self.selected_SPID = diag_ui.getSelectedSPID()
        self.selected_services = diag_ui.getSelectedServices()
        self.selected_tmtc = diag_ui.getTMTC()
        self.showMessage('Applying filter...')

        self.addPacketsToView(self.data, True, show_stat=False)

    def onExportButtonClicked(self):
        if self.y:
            filename = str(
                QtWidgets.QFileDialog.getSaveFileName(
                    None, "Save data to file", "", "CSV(*.csv)")[0])
            if filename:
                with open(filename, 'w') as f:
                    f.write('{},{}\n'.format(self.xlabel, self.ylabel))
                    for xx, yy in zip(self.x, self.y):
                        f.write('{},{}\n'.format(xx, yy))
                    self.showMessage(
                        'The data has been written to {}'.format(filename))
        else:
            msgBox = QtWidgets.QMessageBox()
            msgBox.setIcon(QtWidgets.QMessageBox.Information)
            msgBox.setText('Plot first!')
            msgBox.setWindowTitle("Warning")
            msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
            msgBox.exec_()

    def savePlot(self):
        # if self.figure.get_axes():
        if self.chart:
            filetypes = "PNG (*.png);;JPEG (*.jpg)"
            filename = str(
                QtWidgets.QFileDialog.getSaveFileName(
                    None, "Save plot to file", "", filetypes)[0])
            if filename:
                if not filename.endswith(('.png', '.jpg')):
                    filename += '.png'
                # self.figure.savefig(filename)
                p = self.chartView.grab()
                p.save(filename)
                self.showMessage(('Has been saved to %s.' % filename))

        else:
            msgBox = QtWidgets.QMessageBox()
            msgBox.setIcon(QtWidgets.QMessageBox.Information)
            msgBox.setText('No figure to save')
            msgBox.setWindowTitle("STIX raw data viewer")
            msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
            msgBox.exec_()

    def onCopyTriggered(self):
        packet_id = self.current_row
        try:
            packet = self.data[packet_id]
            ss = pprint.pformat(packet)

            cb = QtWidgets.QApplication.clipboard()
            cb.clear(mode=cb.Clipboard)
            cb.setText(ss, mode=cb.Clipboard)
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Information)
            msg.setText(
                "The data of the selected packet has been copied to the clipboard."
            )
            msg.setWindowTitle("Information")
            msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
            retval = msg.exec_()

        except Exception as e:
            self.showMessage(str(e), 0)

    def onPasteTriggered(self):
        raw_hex = QtWidgets.QApplication.clipboard().text()
        if len(raw_hex) < 16:
            self.showMessage('No data in the clipboard.')
            return
        self.hexParser.setHex(raw_hex)
        self.hexParser.start()

    def showMessage(self, msg, where=0):
        if where != 1:
            self.statusbar.showMessage(msg)
        if where != 0:
            self.statusListWidget.addItem(msg)

    def onSetIDBClicked(self):
        self.idb_filename = QtWidgets.QFileDialog.getOpenFileName(
            None, 'Select file', '.', 'IDB file(*.db *.sqlite *.sqlite3)')[0]

        if not self.idb_filename:
            return

        STIX_IDB.reload(self.idb_filename)
        if STIX_IDB.is_connected():
            #settings = QtCore.QSettings('FHNW', 'parser')
            self.settings.setValue('idb_filename', self.idb_filename)
        self.showMessage(
            'IDB location: {} '.format(STIX_IDB.get_idb_filename()), 1)

    def save(self):
        filetypes = 'python compressed pickle (*.pklz);; python pickle file (*.pkl);; binary data (*.dat)'
        self.output_filename = str(
            QtWidgets.QFileDialog.getSaveFileName(None, "Save packets to", "",
                                                  filetypes)[0])

        if not self.output_filename.endswith(('.pklz', '.pkl', '.dat')):
            msg = 'unsupported file format !'
            self.showMessage(msg)
            return

        msg = 'Writing data to file %s' % self.output_filename
        self.showMessage(msg)
        if self.output_filename.endswith(('.pklz', '.pkl')):
            stw = stixio.StixPickleWriter(self.output_filename)
            stw.register_run(str(self.input_filename))
            stw.write_all(self.data)
        elif self.output_filename.endswith('.dat'):
            stw = stixio.StixBinaryWriter(self.output_filename)
            stw.write_all(self.data)
            num_ok = stw.get_num_sucess()
            msg = (
                'The binary data of {} packets written to file {}, total packets {}'
                .format(num_ok, self.output_filename, len(self.data)))
            self.showMessage(msg)
        msg = 'Packets have been written to %s' % self.output_filename
        self.showMessage(msg)

    def setListViewSelected(self, row):
        #index = self.model.createIndex(row, 0);
        # if index.isValid():
        #    self.model.selectionModel().select( index, QtGui.QItemSelectionModel.Select)
        pass

    def about(self):
        msgBox = QtWidgets.QMessageBox()
        msgBox.setIcon(QtWidgets.QMessageBox.Information)
        msgBox.setText("STIX raw data parser and viewer, hualin.xiao@fhnw.ch")
        msgBox.setWindowTitle("Stix data viewer")

        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msgBox.exec_()

    def nextPacket(self):
        self.current_row += 1
        length = len(self.data)

        if self.current_row >= length:
            self.current_row = length - 1
            self.showMessage('No more packet!')

        self.showPacket(self.current_row)
        self.setListViewSelected(self.current_row)

    def previousPacket(self):
        self.current_row -= 1
        if self.current_row < 0:
            self.current_row = 0
            self.showMessage('Reach the first packet!')

        self.showPacket(self.current_row)
        self.setListViewSelected(self.current_row)

    def getOpenFilename(self):

        location = self.settings.value('location', [], str) 
        if not location:
            location = '.'


        filetypes = (
            'Supported file (*.dat *.bin *.binary *.pkl *.pklz *.xml *ascii *BDF *txt) ;; All(*)'
        )
        self.input_filename = QtWidgets.QFileDialog.getOpenFileName(
            None, 'Select file', location, filetypes)[0]
        if not self.input_filename:
            return
        self.settings.setValue('location', os.path.abspath(self.input_filename))

        diag = QtWidgets.QDialog()
        diag_ui = packet_filter.Ui_Dialog()
        diag_ui.setupUi(diag)
        diag_ui.setSelectedServices(SELECTED_SERVICES)
        diag_ui.buttonBox.accepted.connect(
            partial(self.onOpenFile, self.input_filename, diag_ui))
        diag.exec_()

    def onOpenFile(self, input_filename, diag):
        self.selected_SPID = diag.getSelectedSPID()
        self.selected_services = diag.getSelectedServices()
        self.openFile(input_filename, self.selected_services,
                      self.selected_SPID)

    def openFile(self, filename, selected_services=None, selected_SPID=None):
        msg = 'Loading file %s ...' % filename
        self.progressDiag = QtWidgets.QProgressDialog()
        #self.showMessage(msg)
        self.progressDiag.setLabelText(msg)
        self.progressDiag.setWindowTitle('Loading data')
        self.progressDiag.setCancelButtonText('Cancel')
        self.progressDiag.setRange(0, 100)
        self.progressDiag.setMinimumWidth(300)
        self.progressDiag.canceled.connect(self.stopParsing)
        self.filterPattern = ''
        self.dataReader.setPacketFilter(selected_services, selected_SPID)
        self.dataReader.setFilename(filename)
        self.dataReader.start()
        self.progressDiag.show()

    def stopParsing(self):
        if self.dataReader:
            self.dataReader.stopParsing()
            self.progressDiag.hide()

    def onProgressUpdated(self, progress):
        progress=int(progress)
        if not self.progressDiag:
            return
        self.progressDiag.setValue(progress)
        if progress >=99:
            self.progressDiag.hide()

    def onDataReaderCritical(self, msg):
        self.showMessage(msg, 1)

    def onDataReaderInfo(self, msg):
        self.showMessage(msg, 0)

    def onDataReaderWarning(self, msg):
        self.showMessage(msg, 1)

    def onDataReaderError(self, msg):
        self.showMessage(msg, 1)

    def onDataReady(self, data, clear=True, show_stat=True):
        #self.progressBar.hide()
        if not clear:
            self.data.extend(data)
        else:
            self.data = data
        if data:
            self.addPacketsToView(data, clear=clear, show_stat=show_stat)
            self.enableButtons()
        else:
            self.showMessage('No packet loaded')

    def enableButtons(self):
        if not self.buttons_enabled:
            self.actionPrevious.setEnabled(True)
            self.actionNext.setEnabled(True)
            self.actionSave.setEnabled(True)
            self.actionCopy.setEnabled(True)
            self.actionPlot.setEnabled(True)
            self.actionViewBinary.setEnabled(True)
            self.buttons_enabled = True

    def addPacketsToView(self, data, clear=True, show_stat=True):
        if clear:
            self.packetTreeWidget.clear()

        for p in data:
            if not isinstance(p, dict):
                continue
            header = p['header']
            root = QtWidgets.QTreeWidgetItem(self.packetTreeWidget)
            colors = {2: '#FFA500', 1: '#000080', 3: '#FF0000', 4: '#800000'}
            tc_color = '#78281F'
            if header['TMTC'] == 'TC':
                root.setForeground(0, QtGui.QBrush(QtGui.QColor(tc_color)))
                root.setForeground(1, QtGui.QBrush(QtGui.QColor(tc_color)))
            else:
                if header['service_type'] == 5:
                    if header['service_subtype'] in colors.keys():
                        root.setForeground(
                            0,
                            QtGui.QBrush(
                                QtGui.QColor(
                                    colors[header['service_subtype']])))
                        root.setForeground(
                            1,
                            QtGui.QBrush(
                                QtGui.QColor(
                                    colors[header['service_subtype']])))

            timestamp_str = tu.format_datetime(header['unix_time'])
            root.setText(0, timestamp_str)
            description = '{}({},{}) - {}'.format(
                header['TMTC'], header['service_type'],
                header['service_subtype'], header['descr'])
            root.setText(1, description)
            hidden = False

            if self.selected_SPID:
                if header['TMTC'] == 'TC':
                    hidden = True
                elif -int(header['SPID']) in self.selected_SPID or int(
                        header['SPID']) not in self.selected_SPID:
                    hidden = True
            else:
                if int(header['service_type']) not in self.selected_services:
                    hidden = True

            TMTC = header['TMTC']
            if TMTC == 'TM' and self.selected_tmtc in [2, 0]:
                hidden = True
            if TMTC == 'TC' and self.selected_tmtc in [1, 0]:
                hidden = True


            if self.filterPattern:
                to_exclude = False
                pattern = self.filterPattern.strip()
                if pattern.startswith('!'):
                    to_exclude = True
                    pattern = pattern[1:]
                try:
                    spid = int(pattern)
                    hidden = to_exclude == (header['SPID'] == spid)
                    #XNOR operation
                except (TypeError, ValueError):
                    hidden = to_exclude == (pattern in description)

            root.setHidden(hidden)

        if show_stat:
            total_packets = len(self.data)
            self.showMessage(('Total packet(s): %d' % total_packets))

    def onConnectTSCTriggered(self):
        diag = QtWidgets.QDialog()
        diag_ui = tsc_connection.Ui_Dialog()
        diag_ui.setupUi(diag)
        self.tsc_host = self.settings.value('tsc_host', [], str)
        self.tsc_port = self.settings.value('tsc_port', [], str)
        if self.tsc_host:
            diag_ui.serverLineEdit.setText(self.tsc_host)
        if self.tsc_port:
            diag_ui.portLineEdit.setText(self.tsc_port)
        diag_ui.buttonBox.accepted.connect(partial(self.connectToTSC, diag_ui))
        diag.exec_()

    def connectToTSC(self, dui):
        host = dui.serverLineEdit.text()
        port = dui.portLineEdit.text()
        self.showMessage('Connecting to TSC...')

        self.socketPacketReceiver.connect(host, int(port))
        self.socketPacketReceiver.start()

    def onPacketArrival(self, packets):
        clear = False
        if packets:
            if len(self.data) > MAX_NUM_PACKET_IN_BUFFER:
                clear = True
            self.onDataReady(packets, clear=clear, show_stat=True)

    def onLoadMongoDBTriggered(self):
        diag = QtWidgets.QDialog()
        diag_ui = mongo_dialog.Ui_Dialog()
        diag_ui.setupUi(diag)
        #self.settings = QtCore.QSettings('FHNW', 'parser')
        self.mongo_server = self.settings.value('mongo_server', [], str)
        self.mongo_port = self.settings.value('mongo_port', [], str)
        self.mongo_user = self.settings.value('mongo_user', [], str)
        self.mongo_pwd = self.settings.value('mongo_pwd', [], str)

        if self.mongo_server:
            diag_ui.serverLineEdit.setText(self.mongo_server)
        if self.mongo_port:
            diag_ui.portLineEdit.setText(self.mongo_port)

        if self.mongo_user:
            diag_ui.userLineEdit.setText(self.mongo_user)
        if self.mongo_pwd:
            diag_ui.pwdLineEdit.setText(self.mongo_pwd)

        diag_ui.pushButton.clicked.connect(
            partial(self.loadRunsFromMongoDB, diag_ui))
        diag_ui.buttonBox.accepted.connect(
            partial(self.loadDataFromMongoDB, diag_ui, diag))
        diag.exec_()

    def loadRunsFromMongoDB(self, dui):
        server = dui.serverLineEdit.text()
        port = dui.portLineEdit.text()
        user = dui.userLineEdit.text()
        pwd = dui.pwdLineEdit.text()

        self.showMessage('saving setting...')
        if self.mongo_server != server:
            self.settings.setValue('mongo_server', server)
        if self.mongo_port != port:
            self.settings.setValue('mongo_port', port)

        if self.mongo_user != user:
            self.settings.setValue('mongo_user', user)
        if self.mongo_pwd != pwd:
            self.settings.setValue('mongo_pwd', pwd)

        self.showMessage('connecting Mongo database ...')
        self.mdb = mgdb.MongoDB(server, int(port), user, pwd)
        if not self.mdb.is_connected():
            self.showMessage('Failed to connect to MongoDB')
            return

        dui.treeWidget.clear()
        self.showMessage('Fetching data...')
        for run in self.mdb.select_all_runs():
            root = QtWidgets.QTreeWidgetItem(dui.treeWidget)
            root.setText(0, str(run['_id']))
            root.setText(1, run['filename'])
            root.setText(2, tu.format_datetime(run['date']))
            root.setText(3, tu.format_datetime(run['data_start_unix_time']))
            root.setText(4, tu.format_datetime(run['data_stop_unix_time']))

    def loadDataFromMongoDB(self, dui, diag):
        self.showMessage('Loading packets ...')
        selected_runs = []
        for item in dui.treeWidget.selectedItems():
            selected_runs.append(item.text(0))
        if not selected_runs:
            self.showMessage('Run not selected!')
        if selected_runs:
            diag.done(0)
            self.showMessage('Loading data ...!')
            data = list(self.mdb.select_packets_by_run(selected_runs[0]))
            if data:
                self.onDataReady(data, clear=True)
            else:
                self.showMessage('No packets found!')
        # close

    def onPacketSelected(self, cur, pre):
        self.current_row = self.packetTreeWidget.currentIndex().row()
        self.showMessage(('Packet #%d selected' % self.current_row))
        self.showPacket(self.current_row)

    def showPacket(self, row):
        if not self.data:
            return
        header = self.data[row]['header']
        total_packets = len(self.data)
        self.showMessage(
            ('Packet %d / %d  %s ' % (row, total_packets, header['descr'])))
        self.paramTreeWidget.clear()
        header_root = QtWidgets.QTreeWidgetItem(self.paramTreeWidget)
        header_root.setText(0, "Header")
        rows = len(header)
        for key, val in header.items():
            root = QtWidgets.QTreeWidgetItem(header_root)
            root.setText(0, key)
            root.setText(1, str(val))

        params = self.data[row]['parameters']
        param_root = QtWidgets.QTreeWidgetItem(self.paramTreeWidget)
        param_root.setText(0, "Parameters")
        self.showParameterTree(params, param_root)
        self.paramTreeWidget.expandItem(param_root)
        self.paramTreeWidget.expandItem(header_root)
        self.current_row = row

    def showParameterTree(self, params, parent, parent_id=[]):
        if not params:
            return
        for i, p in enumerate(params):
            root = QtWidgets.QTreeWidgetItem(parent)
            if not p:
                continue
            param = Parameter(p)
            param_name = param['name']
            desc = param['desc']
            current_ids = parent_id[:]
            current_ids.append(i)
            root.setText(0, param_name)

            root.setText(1, desc)
            root.setText(2, str(param['raw']))
            tip='parameter'+''.join(['[{}]'.format(x) for x in current_ids])
            root.setToolTip(0, tip)

            long_desc = STIX_IDB.get_scos_description(param_name)
            if long_desc:
                root.setToolTip(1, long_desc)

            try:
                root.setToolTip(2, hex(param['raw_int']))
            except:
                pass
            unit=STIX_IDB.get_parameter_unit(param_name)
            try:
                eng=int(param['eng'])
            except (TypeError, ValueError):
                eng=''
            try:
                error=int(param['error'])
                eng=f'{eng} \u00B1 {error}'
            except (TypeError, ValueError):
                pass


            root.setText(3, eng)
            root.setText(4, unit)
            if 'NIXG' in param_name:
                root.setHidden(True)
                #groups should not be shown
            if param.children:
                self.showParameterTree(param['children'], root, current_ids)
        self.paramTreeWidget.itemDoubleClicked.connect(self.onTreeItemClicked)

    def walk(self, name, params, header, ret_x, ret_y, xaxis=0, data_type=0):
        if not params:
            return
        timestamp = header['unix_time']
        for p in params:
            if not p:
                continue
            param = Parameter(p)
            if name == param.name:
                values = None
                if data_type == 0:
                    values = param.raw
                else:
                    values = param.eng
                try:
                    yvalue = float(values)
                    ret_y.append(yvalue)
                    if xaxis == 1:
                        ret_x.append(timestamp)
                    else:
                        self.showMessage(('Can not plot %s  ' % str(yvalue)))
                except Exception as e:
                    self.showMessage(('%s ' % str(e)))

            if param.children:
                self.walk(name, param.children, header, ret_x, ret_y, xaxis,
                          data_type)

    def onPlotButtonClicked(self, packets=None):
        if self.chart:
            self.chart.removeAllSeries()
        if packets is None:
            packets = self.data
        if not packets:
            return

        self.showMessage('Preparing plot ...')
        name = self.paramNameEdit.text()
        packet_selection = self.comboBox.currentIndex()
        xaxis_type = self.xaxisComboBox.currentIndex()
        data_type = self.dataTypeComboBox.currentIndex()
        timestamp = []
        self.y = []
        params = self.paramNameEdit.text()
        header = packets[0]['header']
        current_spid = 0
        spid_text = self.spidLineEdit.text()
        if spid_text:
            current_spid = int(spid_text)
        selected_packets=[] 
        if packet_selection == 0:
            selected_packets=[packets[self.current_row]]
        elif packet_selection == 1:
            selected_packets=packets

        for packet in selected_packets:
            header = packet['header']
            if packet['header']['SPID'] != current_spid:
                continue
            params = packet['parameters']
            self.walk(name, params, header, timestamp, self.y, xaxis_type,
                      data_type)


        self.x = []

        if not self.y:
            self.showMessage('No data points')
        elif self.y:
            style = self.styleEdit.text()
            if not style:
                style = '-'
            title = '%s' % str(name)
            desc = self.descLabel.text()
            if desc:
                title += '- %s' % desc

            self.chart.setTitle(title)

            ylabel = 'Raw value'
            xlabel = name
            if data_type == 1:
                ylabel = 'Engineering / Decompressed  value'
            if xaxis_type == 0:
                if packet_selection == 1:
                    xlabel = "Packet #"
                else:
                    xlabel = "Repeat #"
                self.x = range(0, len(self.y))
            if xaxis_type == 1:
                self.x = [t - timestamp[0] for t in timestamp]
                xlabel = 'Time -T0 (s)'

            if xaxis_type != 2:
                series = QLineSeries()
                series2 = None
                for xx, yy in zip(self.x, self.y):
                    series.append(xx, yy)
                if 'o' in style:
                    series2 = QScatterSeries()
                    for xx, yy in zip(self.x, self.y):
                        series2.append(xx, yy)
                    self.chart.addSeries(series2)
                self.chart.addSeries(series)
                axisX = QValueAxis()
                axisX.setTitleText(xlabel)
                axisY = QValueAxis()
                axisY.setTitleText(ylabel)

                self.chart.setAxisX(axisX)
                self.chart.setAxisY(axisY)
                series.attachAxis(axisX)
                series.attachAxis(axisY)
            else:
                nbins = len(set(self.y))
                ycounts, xedges = np.histogram(self.y, bins=nbins)
                series = QLineSeries()
                for i in range(0, nbins):
                    meanx = (xedges[i] + xedges[i + 1]) / 2.
                    series.append(meanx, ycounts[i])
                self.chart.addSeries(series)
                axisX = QValueAxis()
                axisX.setTitleText(name)
                axisY = QValueAxis()
                axisY.setTitleText("Counts")

                self.chart.setAxisY(axisY)
                self.chart.setAxisX(axisX)
                series.attachAxis(axisX)
                series.attachAxis(axisY)

            self.xlabel = xlabel
            self.ylabel = ylabel
            self.chartView.setRubberBand(QChartView.RectangleRubberBand)
            self.chartView.setRenderHint(QtGui.QPainter.Antialiasing)
            msg = 'Number of data points: {}, Ymin: {}, Ymax: {}'.format(
                len(self.y), min(self.y), max(self.y))
            self.showMessage(msg, 1)

            self.showMessage('The canvas updated!')

    def plotParameter(self, SPID=None, pname=None, desc=None):
        self.tabWidget.setCurrentIndex(1)
        if pname:
            self.paramNameEdit.setText(pname)
        if desc:
            self.descLabel.setText(desc)
        if SPID:
            self.spidLineEdit.setText(str(SPID))

    def onPlotActionClicked(self):
        self.tabWidget.setCurrentIndex(1)
        self.plotParameter()

    def onTreeItemClicked(self, it, col):
        SPID = None
        try:
            SPID = self.data[self.current_row]['header']['SPID']
        except IndexError:
            pass
        self.plotParameter(SPID, it.text(0), it.text(1))

