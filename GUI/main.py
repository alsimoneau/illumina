# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main.ui'
#
# Created by: PyQt5 UI code generator 5.15.5
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
from subprocess import call
import yaml
import sys, os, shutil
#import illum
#global variables
switch1 = 0
perc = []
tech = []
ulor = []
inventory_line = []
atm_type = ""

class Ui_ILLUMINA(object):
    def setupUi(self, ILLUMINA):
        ILLUMINA.setObjectName("ILLUMINA")
        ILLUMINA.setEnabled(True)
        ILLUMINA.resize(764, 889)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(ILLUMINA.sizePolicy().hasHeightForWidth())
        ILLUMINA.setSizePolicy(sizePolicy)
        ILLUMINA.setMinimumSize(QtCore.QSize(764, 889))
        ILLUMINA.setMaximumSize(QtCore.QSize(764, 889))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        ILLUMINA.setFont(font)
        self.groupBox_2 = QtWidgets.QGroupBox(ILLUMINA)
        self.groupBox_2.setGeometry(QtCore.QRect(20, 30, 651, 261))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        self.groupBox_2.setMinimumSize(QtCore.QSize(651, 150))
        self.groupBox_2.setMaximumSize(QtCore.QSize(651, 1000))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.groupBox_2.setFont(font)
        self.groupBox_2.setMouseTracking(False)
        self.groupBox_2.setObjectName("groupBox_2")
        self.box_band = QtWidgets.QLineEdit(self.groupBox_2)
        self.box_band.setGeometry(QtCore.QRect(280, 60, 71, 27))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.box_band.setFont(font)
        self.box_band.setReadOnly(True)
        self.box_band.setObjectName("box_band")
        self.box_direction = QtWidgets.QLineEdit(self.groupBox_2)
        self.box_direction.setGeometry(QtCore.QRect(440, 30, 71, 27))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.box_direction.setFont(font)
        self.box_direction.setReadOnly(True)
        self.box_direction.setObjectName("box_direction")
        self.label = QtWidgets.QLabel(self.groupBox_2)
        self.label.setGeometry(QtCore.QRect(10, 32, 80, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setMinimumSize(QtCore.QSize(0, 0))
        self.label.setMaximumSize(QtCore.QSize(80, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label.setFont(font)
        self.label.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Canada))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.groupBox_2)
        self.label_2.setGeometry(QtCore.QRect(10, 62, 85, 21))
        self.label_2.setMaximumSize(QtCore.QSize(85, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.groupBox_2)
        self.label_3.setGeometry(QtCore.QRect(160, 32, 111, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setMaximumSize(QtCore.QSize(111, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.groupBox_2)
        self.label_4.setGeometry(QtCore.QRect(160, 60, 111, 21))
        self.label_4.setMaximumSize(QtCore.QSize(111, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.groupBox_2)
        self.label_5.setGeometry(QtCore.QRect(370, 32, 71, 21))
        self.label_5.setMaximumSize(QtCore.QSize(71, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.latitude_edit = QtWidgets.QLineEdit(self.groupBox_2)
        self.latitude_edit.setGeometry(QtCore.QRect(100, 30, 51, 25))
        self.latitude_edit.setObjectName("latitude_edit")
        self.longitude_edit = QtWidgets.QLineEdit(self.groupBox_2)
        self.longitude_edit.setGeometry(QtCore.QRect(100, 60, 51, 25))
        self.longitude_edit.setObjectName("longitude_edit")
        self.date_edit = QtWidgets.QLineEdit(self.groupBox_2)
        self.date_edit.setGeometry(QtCore.QRect(280, 30, 71, 25))
        self.date_edit.setObjectName("date_edit")
        self.comboBox_atm = QtWidgets.QComboBox(self.groupBox_2)
        self.comboBox_atm.setGeometry(QtCore.QRect(440, 60, 161, 25))
        self.comboBox_atm.setObjectName("comboBox_atm")
        self.comboBox_atm.addItem("")
        self.comboBox_atm.addItem("")
        self.comboBox_atm.addItem("")
        self.comboBox_atm.addItem("")
        self.comboBox_atm.addItem("")
        self.comboBox_atm.addItem("")
        self.comboBox_atm.addItem("")
        self.label_13 = QtWidgets.QLabel(self.groupBox_2)
        self.label_13.setGeometry(QtCore.QRect(370, 60, 71, 21))
        self.label_13.setMaximumSize(QtCore.QSize(71, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.label_11 = QtWidgets.QLabel(self.groupBox_2)
        self.label_11.setGeometry(QtCore.QRect(290, 120, 71, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy)
        self.label_11.setMaximumSize(QtCore.QSize(111, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.comboBox_type = QtWidgets.QComboBox(self.groupBox_2)
        self.comboBox_type.setGeometry(QtCore.QRect(290, 150, 86, 25))
        self.comboBox_type.setObjectName("comboBox_type")
        self.comboBox_type.addItem("")
        self.comboBox_type.addItem("")
        self.comboBox_type.addItem("")
        self.comboBox_type.addItem("")
        self.perc_edit = QtWidgets.QLineEdit(self.groupBox_2)
        self.perc_edit.setGeometry(QtCore.QRect(10, 150, 61, 25))
        self.perc_edit.setObjectName("perc_edit")
        self.comboBox_tech = QtWidgets.QComboBox(self.groupBox_2)
        self.comboBox_tech.setGeometry(QtCore.QRect(90, 150, 86, 25))
        self.comboBox_tech.setObjectName("comboBox_tech")
        self.comboBox_tech.addItem("")
        self.comboBox_tech.addItem("")
        self.comboBox_tech.addItem("")
        self.comboBox_tech.addItem("")
        self.comboBox_tech.addItem("")
        self.comboBox_tech.addItem("")
        self.comboBox_tech.addItem("")
        self.label_10 = QtWidgets.QLabel(self.groupBox_2)
        self.label_10.setGeometry(QtCore.QRect(10, 120, 111, 21))
        self.label_10.setMaximumSize(QtCore.QSize(111, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.comboBox_ulor = QtWidgets.QComboBox(self.groupBox_2)
        self.comboBox_ulor.setGeometry(QtCore.QRect(190, 150, 86, 25))
        self.comboBox_ulor.setObjectName("comboBox_ulor")
        self.comboBox_ulor.addItem("")
        self.comboBox_ulor.addItem("")
        self.comboBox_ulor.addItem("")
        self.comboBox_ulor.addItem("")
        self.comboBox_ulor.addItem("")
        self.comboBox_ulor.addItem("")
        self.comboBox_ulor.addItem("")
        self.label_12 = QtWidgets.QLabel(self.groupBox_2)
        self.label_12.setGeometry(QtCore.QRect(200, 120, 61, 21))
        self.label_12.setMaximumSize(QtCore.QSize(111, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.label_9 = QtWidgets.QLabel(self.groupBox_2)
        self.label_9.setGeometry(QtCore.QRect(90, 120, 111, 21))
        self.label_9.setMaximumSize(QtCore.QSize(111, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.add_tech_btn = QtWidgets.QPushButton(self.groupBox_2)
        self.add_tech_btn.setGeometry(QtCore.QRect(400, 150, 61, 25))
        self.add_tech_btn.setObjectName("add_tech_btn")
        self.add_source_btn = QtWidgets.QPushButton(self.groupBox_2)
        self.add_source_btn.setGeometry(QtCore.QRect(470, 150, 141, 25))
        self.add_source_btn.setObjectName("add_source_btn")
        self.exp_def_button = QtWidgets.QPushButton(self.groupBox_2)
        self.exp_def_button.setGeometry(QtCore.QRect(260, 210, 121, 25))
        self.exp_def_button.setObjectName("exp_def_button")

        self.retranslateUi(ILLUMINA)
        QtCore.QMetaObject.connectSlotsByName(ILLUMINA)

        self.exp_def_button.clicked.connect(self.defining_exp)
        self.add_source_btn.clicked.connect(self.defining_source)
        self.add_tech_btn.clicked.connect(self.defining_combo)

    def retranslateUi(self, ILLUMINA):
        _translate = QtCore.QCoreApplication.translate
        ILLUMINA.setWindowTitle(_translate("ILLUMINA", "ILLUMINA LIGHT"))
        self.groupBox_2.setTitle(_translate("ILLUMINA", "Input definition"))
        self.box_band.setText(_translate("ILLUMINA", "V (J-C)"))
        self.box_direction.setText(_translate("ILLUMINA", "Zenith"))
        self.label.setText(_translate("ILLUMINA", "Latitude (deg)"))
        self.label_2.setText(_translate("ILLUMINA", "Longitude (deg)"))
        self.label_3.setText(_translate("ILLUMINA", "Date (dd/mm/yyyy)"))
        self.label_4.setText(_translate("ILLUMINA", "Band"))
        self.label_5.setText(_translate("ILLUMINA", "Direction"))
        self.comboBox_atm.setItemText(0, _translate("ILLUMINA", "CC Continental clean"))
        self.comboBox_atm.setItemText(1, _translate("ILLUMINA", "CA Continental average"))
        self.comboBox_atm.setItemText(2, _translate("ILLUMINA", "CP Continental polluted"))
        self.comboBox_atm.setItemText(3, _translate("ILLUMINA", "MP Maritime polluted"))
        self.comboBox_atm.setItemText(4, _translate("ILLUMINA", "MC Maritime clean"))
        self.comboBox_atm.setItemText(5, _translate("ILLUMINA", "U Urban"))
        self.comboBox_atm.setItemText(6, _translate("ILLUMINA", "D Desert"))
        self.label_13.setText(_translate("ILLUMINA", "Atmosphere"))
        self.label_11.setText(_translate("ILLUMINA", "Type"))
        self.comboBox_type.setCurrentText(_translate("ILLUMINA", "City"))
        self.comboBox_type.setItemText(0, _translate("ILLUMINA", "City"))
        self.comboBox_type.setItemText(1, _translate("ILLUMINA", "Town"))
        self.comboBox_type.setItemText(2, _translate("ILLUMINA", "Village"))
        self.comboBox_type.setItemText(3, _translate("ILLUMINA", "Rural"))
        self.comboBox_tech.setItemText(0, _translate("ILLUMINA", "HPS"))
        self.comboBox_tech.setItemText(1, _translate("ILLUMINA", "3LED"))
        self.comboBox_tech.setItemText(2, _translate("ILLUMINA", "4LED"))
        self.comboBox_tech.setItemText(3, _translate("ILLUMINA", "LED PC Ambre"))
        self.comboBox_tech.setItemText(4, _translate("ILLUMINA", "MV"))
        self.comboBox_tech.setItemText(5, _translate("ILLUMINA", "CFL"))
        self.comboBox_tech.setItemText(6, _translate("ILLUMINA", "MH"))
        self.label_10.setText(_translate("ILLUMINA", "%"))
        self.comboBox_ulor.setItemText(0, _translate("ILLUMINA", "0"))
        self.comboBox_ulor.setItemText(1, _translate("ILLUMINA", "1"))
        self.comboBox_ulor.setItemText(2, _translate("ILLUMINA", "5"))
        self.comboBox_ulor.setItemText(3, _translate("ILLUMINA", "10"))
        self.comboBox_ulor.setItemText(4, _translate("ILLUMINA", "15"))
        self.comboBox_ulor.setItemText(5, _translate("ILLUMINA", "20"))
        self.comboBox_ulor.setItemText(6, _translate("ILLUMINA", "50"))
        self.label_12.setText(_translate("ILLUMINA", "ULOR"))
        self.label_9.setText(_translate("ILLUMINA", "Technology"))
        self.add_tech_btn.setText(_translate("ILLUMINA", "Add"))
        self.add_source_btn.setText(_translate("ILLUMINA", "Finished"))
        self.exp_def_button.setText(_translate("ILLUMINA", "Create Inputs"))

    def defining_exp(self):
        global atm_type
        latitude = self.latitude_edit.text()
        longitude = self.longitude_edit.text()
        # date = self.date_edit.text()
        # date_day = date.split('/')[0]
        # date_month = date.split('/')[1]
        # date_year = date.split('/')[2]

        self.box_band.setText('V J-C')
        self.box_direction.setText('Zenith')
        atm_type = self.comboBox_atm.currentText()
        with open('domain_params.in', 'w') as f:
            f.write('latitude: ' + latitude + '\n')
            f.write('longitude: ' + longitude + '\n')
            f.write('srs: auto\n')
            f.write('scale_factor: 3\n')
            f.write('nb_pixels: 10\n')
            f.write('nb_layers: 3\n')
            f.write('scale_min: 1000\n')
            f.write('buffer: 10')
        global inventory_line
        if os.path.isfile('inventory.txt'):
            os.remove('inventory.txt')
        with open('inventory.txt', 'w') as f:
            for line in inventory_line:
                f.write(line)
        inventory_line = []
        call(["illum","domain"])
        call(["illum","warp"])
        with open('inputs_params.in', 'r') as f:
            lines = f.readlines()

        with open('inputs_params.in', 'w') as f:
            for line in lines:
                if line[:15] == 'aerosol_profile':
                    f.write('aerosol_profile: ' + atm_type.split(' ')[0] + '\n')
                else:
                    f.write(line)

        call(["illum","inputs"])
        source = "./Inputs"
        os.chdir(source)
        #destination = "./"
        #os.remove('./Inputs/srtm.hdf5')
        #files = os.listdir(source)
        #for file in files:
            #shutil.move(f"{source}/{file}", destination)
        call(["illum","batches"])

    def defining_source(self):
        global switch1
        global perc
        global tech
        global ulor
        global inventory_line
        switch1 = 0
        lat = self.latitude_edit.text()
        long = self.longitude_edit.text()
        radius = str(70)
        type = self.comboBox_type.currentText()
        if type == 'City':
            hobs = str(0)
            dobs = str(6)
            fobs = str(0)
            hlamp = str(6)
        elif type == 'Rural':
            hobs = str(5)
            dobs = str(4)
            fobs = str(0.5)
            hlamp = str(5)

        a1 = lat + '\t' + long + '\t' + radius + '\t' + hobs + '\t' + dobs + '\t' + fobs + '\t' + hlamp + '\t'
        a = []
        for i in range(len(perc)):
            if i == len(perc) - 1:
                a.append(perc[i] + '_' + tech[i] + '_' + ulor[i] + '\n')
            else:
                a.append(perc[i] + '_' + tech[i] + '_' + ulor[i] + ' ')
        a2 = ""
        for combo in a:
            a2 += str(combo)
        inventory_line.append(a1 + a2)
        self.perc_edit.setText('')
        self.comboBox_type.setCurrentIndex(0)
        self.comboBox_tech.setCurrentIndex(0)
        self.comboBox_ulor.setCurrentIndex(0)

    def defining_combo(self):
        global switch1
        global perc
        global tech
        global ulor
        if switch1 == 0:
            perc = []
            tech = []
            ulor = []
        perc.append(self.perc_edit.text())
        tech.append(self.comboBox_tech.currentText())
        ulor.append(self.comboBox_ulor.currentText())
        switch1 = 1
        self.perc_edit.setText('')
        self.comboBox_tech.setCurrentIndex(0)
        self.comboBox_ulor.setCurrentIndex(0)


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    ILLUMINA = QtWidgets.QWidget()
    ui = Ui_ILLUMINA()
    ui.setupUi(ILLUMINA)
    ILLUMINA.show()
    sys.exit(app.exec_())
