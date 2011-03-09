# -*- coding: utf-8 -*-

from PyQt4 import QtGui
from PyQt4 import QtCore
from PyQt4.QtCore import SIGNAL, SLOT
from ui_mainwindow import *
#from world import WorldArea

class MainWindow(QtGui.QMainWindow, Ui_MainWindow):

    # Coordinate system mode
    mode = 'cartesian'
    mode = 'spherical'

    # Quantum numbers
    n = 1
    l = 0
    m = 0

    # Cartesian limits and values
    # Lower
    min_x_value = -1
    min_y_value = -1
    min_z_value = -1
    lower_x_limit = -100
    lower_y_limit = -100
    lower_z_limit = -100
    # Upper
    max_x_value = 1
    max_y_value = 1
    max_z_value = 1
    lower_x_limit = 1
    lower_y_limit = 1
    lower_z_limit = 1

    #Spherical limits and values
    min_r = 0

    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self, parent)

        self.setupUi(self)
        self.saveAsActs = []

        self.createActions()
        self.createMenus()

        self.setWindowTitle(self.tr("Molecula"))

        self.connect(self.spin_n, SIGNAL('valueChanged(int)'),
                      self, SLOT('set_n(int)'))
        self.connect(self.spin_l, SIGNAL('valueChanged(int)'),
                      self, SLOT('set_l(int)'))
        self.connect(self.spin_m, SIGNAL('valueChanged(int)'),
                      self, SLOT('set_m(int)'))


        self.connect(self.spinMin_x_r, SIGNAL('valueChanged(int)'),
                      self, SLOT('changedMinimum_x_r(int)'))
        self.connect(self.spinMin_y_theta, SIGNAL('valueChanged(int)'),
                      self, SLOT('changedMinimum_y_theta(int)'))
        self.connect(self.spinMin_z_phi, SIGNAL('valueChanged(int)'),
                      self, SLOT('changedMinimum_z_phi(int)'))

        self.connect(self.cmbSysCoord, SIGNAL('currentIndexChanged(const QString)'), 
                      self, SLOT('toggleMode(const QString)'))

        mode_id = self.cmbSysCoord.findText(self.mode)
        if mode_id > 0:
            self.cmbSysCoord.setCurrentIndex(mode_id)
            self.cmbSysCoord.emit(SIGNAL('currentIndexChanged(const QString)'), self.mode)

        self.spin_n.setValue(self.n)
        self.spin_n.emit(SIGNAL('valueChanged(int)'), self.n)

    @QtCore.pyqtSlot('int')
    def set_n(self, new_value):
        self.n = new_value
        if new_value <= self.l:
            self.l = new_value - 1
        self.spin_l.setMaximum(new_value - 1)
        self.spin_m.setMaximum(self.l)
        self.spin_m.setMinimum(-self.l)

    @QtCore.pyqtSlot('int')
    def set_l(self, new_value):
        self.l = new_value
        self.spin_m.setMaximum(self.l)
        self.spin_m.setMinimum(-self.l)

    @QtCore.pyqtSlot('int')
    def set_m(self, new_value):
        self.m = new_value

    @QtCore.pyqtSlot('int')
    def changedMinimum_x_r(self, new_value):
        self.spinMax_x_r.setMinimum(new_value)

    @QtCore.pyqtSlot('int')
    def changedMinimum_y_theta(self, new_value):
        self.spinMax_y_theta.setMinimum(new_value)

    @QtCore.pyqtSlot('int')
    def changedMinimum_z_phi(self, new_value):
        self.spinMax_z_phi.setMinimum(new_value)

    @QtCore.pyqtSlot('QString')
    def toggleMode(self, mode):
        if mode == self.mode:
            return 0
        else:
            self.mode = mode

        if mode == 'cartesian':
            self.lblMin_x_r.setText('x')
            self.lblMin_x_r.setMinimum(self.min_x)
            self.lblMin_y_theta.setText('y')
            self.lblMin_z_phi.setText('z')
            self.lblMax_x_r.setText('x')
            self.lblMax_y_theta.setText('y')
            self.lblMax_z_phi.setText('z')
        elif mode == 'spherical':
            self.lblMin_x_r.setText('r')
            self.lblMin_x_r.setMinimum(self.min_r)

            self.lblMin_y_theta.setText(u'θ')
            self.lblMin_z_phi.setText(u'φ')

            self.lblMax_x_r.setText('r')
            self.lblMax_x_r.setMinimum(self.min_r)

            self.lblMax_y_theta.setText(u'θ')
            self.lblMax_z_phi.setText(u'φ')
            
    def closeEvent(self, event):
        #if self.maybeSave():
        #    event.accept()
        #else:
        #    event.ignore()
        event.accept()

    def open(self):
        return
    #    if self.maybeSave():
    #        fileName = QtGui.QFileDialog.getOpenFileName(self,
    #                                                     self.tr("Open File"),
    #                                                     QtCore.QDir.currentPath())
    #        if not fileName.isEmpty():
    #            self.world.openImage(fileName)

    def save(self):
        action = self.sender()
        fileFormat = action.data().toByteArray()
        self.saveFile(fileFormat)

    def penColor(self):
        return
    #    newColor = QtGui.QColorDialog.getColor(self.world.penColor())
    #    if newColor.isValid():
    #        self.world.setPenColor(newColor)

    def penWidth(self):
        return
    #    newWidth, ok = QtGui.QInputDialog.getInteger(self, self.tr("Life"),
    #                                           self.tr("Select pen width:"),
    #                                           self.world.penWidth(),
    #                                           1, 50, 1)
    #    if ok:
    #        self.world.setPenWidth(newWidth)

    def about(self):
        QtGui.QMessageBox.about(self, self.tr("About Life"), self.tr(
          "<p>"
          "</p>"))

    def createActions(self):
        self.openAct = QtGui.QAction(self.tr("&Open..."), self)
        self.openAct.setShortcut(self.tr("Ctrl+O"))
        self.connect(self.openAct, QtCore.SIGNAL("triggered()"), self.open)

        for format in QtGui.QImageWriter.supportedImageFormats():
            text = self.tr("%1...").arg(QtCore.QString(format).toUpper())

            action = QtGui.QAction(text, self)
            action.setData(QtCore.QVariant(format))
            self.connect(action, QtCore.SIGNAL("triggered()"), self.save)
            self.saveAsActs.append(action)

        self.exitAct = QtGui.QAction(self.tr("E&xit"), self)
        self.exitAct.setShortcut(self.tr("Ctrl+Q"))
        self.connect(self.exitAct, QtCore.SIGNAL("triggered()"),
                     self, QtCore.SLOT("close()"))

        self.penColorAct = QtGui.QAction(self.tr("&Pen Color..."), self)
        self.connect(self.penColorAct, QtCore.SIGNAL("triggered()"),
                     self.penColor)

        self.penWidthAct = QtGui.QAction(self.tr("Pen &Width..."), self)
        self.connect(self.penWidthAct, QtCore.SIGNAL("triggered()"),
                     self.penWidth)

        #self.clearScreenAct = QtGui.QAction(self.tr("&Clear Screen"), self)
        #self.clearScreenAct.setShortcut(self.tr("Ctrl+L"))
        #self.connect(self.clearScreenAct, QtCore.SIGNAL("triggered()"),
        #             self.world.clearImage)

        self.aboutAct = QtGui.QAction(self.tr("&About"), self)
        self.connect(self.aboutAct, QtCore.SIGNAL("triggered()"), self.about)

        self.aboutQtAct = QtGui.QAction(self.tr("About &Qt"), self)
        self.connect(self.aboutQtAct, QtCore.SIGNAL("triggered()"),
                     QtGui.qApp, QtCore.SLOT("aboutQt()"))

    def createMenus(self):
        self.saveAsMenu = QtGui.QMenu(self.tr("&Save As"), self)
        for action in self.saveAsActs:
            self.saveAsMenu.addAction(action)

        self.fileMenu = QtGui.QMenu(self.tr("&File"), self)
        self.fileMenu.addAction(self.openAct)
        self.fileMenu.addMenu(self.saveAsMenu)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)

        self.optionMenu = QtGui.QMenu(self.tr("&Options"), self)
        self.optionMenu.addAction(self.penColorAct)
        self.optionMenu.addAction(self.penWidthAct)
        self.optionMenu.addSeparator()
        #self.optionMenu.addAction(self.clearScreenAct)

        self.helpMenu = QtGui.QMenu(self.tr("&Help"), self)
        self.helpMenu.addAction(self.aboutAct)
        self.helpMenu.addAction(self.aboutQtAct)

        self.menuBar().addMenu(self.fileMenu)
        self.menuBar().addMenu(self.optionMenu)
        self.menuBar().addMenu(self.helpMenu)

    def maybeSave(self):
        #if self.world.isModified():
        #    ret = QtGui.QMessageBox.warning(self, self.tr("Life"),
        #                self.tr("The image has been modified.\n"
        #                        "Do you want to save your changes?"),
        #                QtGui.QMessageBox.Yes | QtGui.QMessageBox.Default,
        #                QtGui.QMessageBox.No,
        #                QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Escape)
        #    if ret == QtGui.QMessageBox.Yes:
        #        return self.saveFile("png")
        #    elif ret == QtGui.QMessageBox.Cancel:
        #        return False

        return True

    def saveFile(self, fileFormat):
        initialPath = QtCore.QDir.currentPath() + "/untitled." + fileFormat

        fileName = QtGui.QFileDialog.getSaveFileName(self, self.tr("Save As"),
                                    initialPath,
                                    self.tr("%1 Files (*.%2);;All Files (*)")
                                    .arg(QtCore.QString(fileFormat.toUpper()))
                                    .arg(QtCore.QString(fileFormat)))
        if fileName.isEmpty():
            return False
        #else:
        #    return self.world.saveImage(fileName, fileFormat)


