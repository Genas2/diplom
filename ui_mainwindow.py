# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainwindow.ui'
#
# Created: Tue Mar 22 08:36:21 2011
#      by: PyQt4 UI code generator 4.8.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(287, 357)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.verticalLayout_14 = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout_14.setObjectName(_fromUtf8("verticalLayout_14"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.lbl_n = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setItalic(True)
        self.lbl_n.setFont(font)
        self.lbl_n.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_n.setObjectName(_fromUtf8("lbl_n"))
        self.verticalLayout.addWidget(self.lbl_n)
        self.spin_n = QtGui.QSpinBox(self.centralwidget)
        self.spin_n.setMinimum(1)
        self.spin_n.setObjectName(_fromUtf8("spin_n"))
        self.verticalLayout.addWidget(self.spin_n)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.lbl_l = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setPointSize(11)
        font.setItalic(True)
        self.lbl_l.setFont(font)
        self.lbl_l.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_l.setObjectName(_fromUtf8("lbl_l"))
        self.verticalLayout_2.addWidget(self.lbl_l)
        self.spin_l = QtGui.QSpinBox(self.centralwidget)
        self.spin_l.setObjectName(_fromUtf8("spin_l"))
        self.verticalLayout_2.addWidget(self.spin_l)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.lbl_m = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setPointSize(11)
        font.setItalic(True)
        self.lbl_m.setFont(font)
        self.lbl_m.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_m.setObjectName(_fromUtf8("lbl_m"))
        self.verticalLayout_3.addWidget(self.lbl_m)
        self.spin_m = QtGui.QSpinBox(self.centralwidget)
        self.spin_m.setObjectName(_fromUtf8("spin_m"))
        self.verticalLayout_3.addWidget(self.spin_m)
        self.horizontalLayout.addLayout(self.verticalLayout_3)
        self.verticalLayout_14.addLayout(self.horizontalLayout)
        self.verticalLayout_12 = QtGui.QVBoxLayout()
        self.verticalLayout_12.setObjectName(_fromUtf8("verticalLayout_12"))
        self.cmbSysCoord = QtGui.QComboBox(self.centralwidget)
        self.cmbSysCoord.setObjectName(_fromUtf8("cmbSysCoord"))
        self.cmbSysCoord.addItem(_fromUtf8(""))
        self.cmbSysCoord.addItem(_fromUtf8(""))
        self.verticalLayout_12.addWidget(self.cmbSysCoord)
        self.verticalLayout_7 = QtGui.QVBoxLayout()
        self.verticalLayout_7.setObjectName(_fromUtf8("verticalLayout_7"))
        self.lbl_x_r_2 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setItalic(True)
        self.lbl_x_r_2.setFont(font)
        self.lbl_x_r_2.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_x_r_2.setObjectName(_fromUtf8("lbl_x_r_2"))
        self.verticalLayout_7.addWidget(self.lbl_x_r_2)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.lblMin_x_r = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setItalic(True)
        self.lblMin_x_r.setFont(font)
        self.lblMin_x_r.setAlignment(QtCore.Qt.AlignCenter)
        self.lblMin_x_r.setObjectName(_fromUtf8("lblMin_x_r"))
        self.verticalLayout_4.addWidget(self.lblMin_x_r)
        self.spinMin_x_r = QtGui.QDoubleSpinBox(self.centralwidget)
        self.spinMin_x_r.setMinimum(-99.99)
        self.spinMin_x_r.setProperty(_fromUtf8("value"), -1.0)
        self.spinMin_x_r.setObjectName(_fromUtf8("spinMin_x_r"))
        self.verticalLayout_4.addWidget(self.spinMin_x_r)
        self.horizontalLayout_2.addLayout(self.verticalLayout_4)
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setObjectName(_fromUtf8("verticalLayout_5"))
        self.lblMin_y_theta = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setPointSize(11)
        font.setItalic(True)
        self.lblMin_y_theta.setFont(font)
        self.lblMin_y_theta.setAlignment(QtCore.Qt.AlignCenter)
        self.lblMin_y_theta.setObjectName(_fromUtf8("lblMin_y_theta"))
        self.verticalLayout_5.addWidget(self.lblMin_y_theta)
        self.spinMin_y_theta = QtGui.QDoubleSpinBox(self.centralwidget)
        self.spinMin_y_theta.setMinimum(-99.99)
        self.spinMin_y_theta.setProperty(_fromUtf8("value"), -1.0)
        self.spinMin_y_theta.setObjectName(_fromUtf8("spinMin_y_theta"))
        self.verticalLayout_5.addWidget(self.spinMin_y_theta)
        self.horizontalLayout_2.addLayout(self.verticalLayout_5)
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setObjectName(_fromUtf8("verticalLayout_6"))
        self.lblMin_z_phi = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setPointSize(11)
        font.setItalic(True)
        self.lblMin_z_phi.setFont(font)
        self.lblMin_z_phi.setAlignment(QtCore.Qt.AlignCenter)
        self.lblMin_z_phi.setObjectName(_fromUtf8("lblMin_z_phi"))
        self.verticalLayout_6.addWidget(self.lblMin_z_phi)
        self.spinMin_z_phi = QtGui.QDoubleSpinBox(self.centralwidget)
        self.spinMin_z_phi.setMinimum(-99.99)
        self.spinMin_z_phi.setProperty(_fromUtf8("value"), -1.0)
        self.spinMin_z_phi.setObjectName(_fromUtf8("spinMin_z_phi"))
        self.verticalLayout_6.addWidget(self.spinMin_z_phi)
        self.horizontalLayout_2.addLayout(self.verticalLayout_6)
        self.verticalLayout_7.addLayout(self.horizontalLayout_2)
        self.verticalLayout_12.addLayout(self.verticalLayout_7)
        self.verticalLayout_8 = QtGui.QVBoxLayout()
        self.verticalLayout_8.setObjectName(_fromUtf8("verticalLayout_8"))
        self.lbl_x_r_3 = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setItalic(True)
        self.lbl_x_r_3.setFont(font)
        self.lbl_x_r_3.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_x_r_3.setObjectName(_fromUtf8("lbl_x_r_3"))
        self.verticalLayout_8.addWidget(self.lbl_x_r_3)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.verticalLayout_9 = QtGui.QVBoxLayout()
        self.verticalLayout_9.setObjectName(_fromUtf8("verticalLayout_9"))
        self.lblMax_x_r = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setItalic(True)
        self.lblMax_x_r.setFont(font)
        self.lblMax_x_r.setAlignment(QtCore.Qt.AlignCenter)
        self.lblMax_x_r.setObjectName(_fromUtf8("lblMax_x_r"))
        self.verticalLayout_9.addWidget(self.lblMax_x_r)
        self.spinMax_x_r = QtGui.QDoubleSpinBox(self.centralwidget)
        self.spinMax_x_r.setMinimum(-99.99)
        self.spinMax_x_r.setProperty(_fromUtf8("value"), 1.0)
        self.spinMax_x_r.setObjectName(_fromUtf8("spinMax_x_r"))
        self.verticalLayout_9.addWidget(self.spinMax_x_r)
        self.horizontalLayout_4.addLayout(self.verticalLayout_9)
        self.verticalLayout_10 = QtGui.QVBoxLayout()
        self.verticalLayout_10.setObjectName(_fromUtf8("verticalLayout_10"))
        self.lblMax_y_theta = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setPointSize(11)
        font.setItalic(True)
        self.lblMax_y_theta.setFont(font)
        self.lblMax_y_theta.setAlignment(QtCore.Qt.AlignCenter)
        self.lblMax_y_theta.setObjectName(_fromUtf8("lblMax_y_theta"))
        self.verticalLayout_10.addWidget(self.lblMax_y_theta)
        self.spinMax_y_theta = QtGui.QDoubleSpinBox(self.centralwidget)
        self.spinMax_y_theta.setMinimum(-99.99)
        self.spinMax_y_theta.setProperty(_fromUtf8("value"), 1.0)
        self.spinMax_y_theta.setObjectName(_fromUtf8("spinMax_y_theta"))
        self.verticalLayout_10.addWidget(self.spinMax_y_theta)
        self.horizontalLayout_4.addLayout(self.verticalLayout_10)
        self.verticalLayout_11 = QtGui.QVBoxLayout()
        self.verticalLayout_11.setObjectName(_fromUtf8("verticalLayout_11"))
        self.lblMax_z_phi = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Courier New"))
        font.setPointSize(11)
        font.setItalic(True)
        self.lblMax_z_phi.setFont(font)
        self.lblMax_z_phi.setAlignment(QtCore.Qt.AlignCenter)
        self.lblMax_z_phi.setObjectName(_fromUtf8("lblMax_z_phi"))
        self.verticalLayout_11.addWidget(self.lblMax_z_phi)
        self.spinMax_z_phi = QtGui.QDoubleSpinBox(self.centralwidget)
        self.spinMax_z_phi.setMinimum(-99.99)
        self.spinMax_z_phi.setProperty(_fromUtf8("value"), 1.0)
        self.spinMax_z_phi.setObjectName(_fromUtf8("spinMax_z_phi"))
        self.verticalLayout_11.addWidget(self.spinMax_z_phi)
        self.horizontalLayout_4.addLayout(self.verticalLayout_11)
        self.verticalLayout_8.addLayout(self.horizontalLayout_4)
        self.verticalLayout_12.addLayout(self.verticalLayout_8)
        self.verticalLayout_14.addLayout(self.verticalLayout_12)
        self.verticalLayout_13 = QtGui.QVBoxLayout()
        self.verticalLayout_13.setObjectName(_fromUtf8("verticalLayout_13"))
        self.cmbEquations = QtGui.QComboBox(self.centralwidget)
        self.cmbEquations.setObjectName(_fromUtf8("cmbEquations"))
        self.cmbEquations.addItem(_fromUtf8(""))
        self.verticalLayout_13.addWidget(self.cmbEquations)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.btnBuildPlot = QtGui.QPushButton(self.centralwidget)
        self.btnBuildPlot.setObjectName(_fromUtf8("btnBuildPlot"))
        self.horizontalLayout_3.addWidget(self.btnBuildPlot)
        self.btnPreview = QtGui.QPushButton(self.centralwidget)
        self.btnPreview.setObjectName(_fromUtf8("btnPreview"))
        self.horizontalLayout_3.addWidget(self.btnPreview)
        self.btnExit = QtGui.QPushButton(self.centralwidget)
        self.btnExit.setObjectName(_fromUtf8("btnExit"))
        self.horizontalLayout_3.addWidget(self.btnExit)
        self.verticalLayout_13.addLayout(self.horizontalLayout_3)
        self.verticalLayout_14.addLayout(self.verticalLayout_13)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 287, 22))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.btnExit, QtCore.SIGNAL(_fromUtf8("clicked()")), MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.lbl_n.setText(QtGui.QApplication.translate("MainWindow", "n", None, QtGui.QApplication.UnicodeUTF8))
        self.lbl_l.setText(QtGui.QApplication.translate("MainWindow", "l", None, QtGui.QApplication.UnicodeUTF8))
        self.lbl_m.setText(QtGui.QApplication.translate("MainWindow", "m", None, QtGui.QApplication.UnicodeUTF8))
        self.cmbSysCoord.setItemText(0, QtGui.QApplication.translate("MainWindow", "cartesian", None, QtGui.QApplication.UnicodeUTF8))
        self.cmbSysCoord.setItemText(1, QtGui.QApplication.translate("MainWindow", "spherical", None, QtGui.QApplication.UnicodeUTF8))
        self.lbl_x_r_2.setText(QtGui.QApplication.translate("MainWindow", "lower limits", None, QtGui.QApplication.UnicodeUTF8))
        self.lblMin_x_r.setText(QtGui.QApplication.translate("MainWindow", "x", None, QtGui.QApplication.UnicodeUTF8))
        self.lblMin_y_theta.setText(QtGui.QApplication.translate("MainWindow", "y", None, QtGui.QApplication.UnicodeUTF8))
        self.lblMin_z_phi.setText(QtGui.QApplication.translate("MainWindow", "z", None, QtGui.QApplication.UnicodeUTF8))
        self.lbl_x_r_3.setText(QtGui.QApplication.translate("MainWindow", "upper limits", None, QtGui.QApplication.UnicodeUTF8))
        self.lblMax_x_r.setText(QtGui.QApplication.translate("MainWindow", "x", None, QtGui.QApplication.UnicodeUTF8))
        self.lblMax_y_theta.setText(QtGui.QApplication.translate("MainWindow", "y", None, QtGui.QApplication.UnicodeUTF8))
        self.lblMax_z_phi.setText(QtGui.QApplication.translate("MainWindow", "z", None, QtGui.QApplication.UnicodeUTF8))
        self.cmbEquations.setItemText(0, QtGui.QApplication.translate("MainWindow", "Angular part", None, QtGui.QApplication.UnicodeUTF8))
        self.btnBuildPlot.setText(QtGui.QApplication.translate("MainWindow", "Build Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.btnPreview.setText(QtGui.QApplication.translate("MainWindow", "Preview equation", None, QtGui.QApplication.UnicodeUTF8))
        self.btnExit.setText(QtGui.QApplication.translate("MainWindow", "Exit", None, QtGui.QApplication.UnicodeUTF8))

