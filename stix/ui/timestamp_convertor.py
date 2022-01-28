# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'timestamp_convertor.ui'
#
# Created by: PyQt5 UI code generator 5.13.2
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets
from stix.spice import datetime 
from PyQt5.QtCore import QDateTime,Qt

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(466, 232)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.comboBox = QtWidgets.QComboBox(Dialog)
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.gridLayout.addWidget(self.comboBox, 0, 0, 1, 1)
        self.lineEdit = QtWidgets.QLineEdit(Dialog)
        self.lineEdit.setObjectName("lineEdit")
        self.gridLayout.addWidget(self.lineEdit, 0, 1, 1, 1)
        self.pushButton = QtWidgets.QPushButton(Dialog)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 0, 2, 1, 1)
        self.textEdit = QtWidgets.QTextEdit(Dialog)
        self.textEdit.setObjectName("textEdit")
        self.gridLayout.addWidget(self.textEdit, 1, 0, 1, 3)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 2, 1, 1, 2)
        self.verticalLayout.addLayout(self.gridLayout)
        now = QDateTime.currentDateTime()
        self.lineEdit.setText(now.toUTC().toString(Qt.ISODate)[0:-1])
        self.pushButton.clicked.connect(self.convertTimestamp)

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def convertTimestamp(self):
        inputType= str(self.comboBox.currentText())
        inputData=str(self.lineEdit.text())
        outputText=''
        utc=''
        scet=''
        unix=''
        try:
            if inputType == 'UTC':
                utc=inputData
                unix=datetime.utc2unix(inputData)
                scet=datetime.utc2scet(inputData)
            elif inputType=='SCET':
                scet=int(inputData)
                unix=datetime.scet2unix(inputData)
                utc=datetime.scet2utc(inputData)
            else:
                unix=int(inputData)
                utc=datetime.unix2utc(unix)
                scet=datetime.utc2scet(utc)
                
            gap=946728000
            J2000=unix-gap
            scet_utc=datetime.unix2utc(scet+gap)
            outputText='''UTC: {}\nUNIX:{}\nSCET:{}\nJ2000 (calculated from UNIX epoch):{}\n UTC(calculated from SCET):{}
            '''.format(utc, unix,scet,J2000, scet_utc)
        except Exception as e:
            outputText=str(e)

        self.textEdit.setText(outputText)



    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Timestamp convertor"))
        self.comboBox.setItemText(0, _translate("Dialog", "UTC"))
        self.comboBox.setItemText(1, _translate("Dialog", "SCET"))
        self.comboBox.setItemText(2, _translate("Dialog", "Unix Time"))
        self.pushButton.setText(_translate("Dialog", "Convert"))
