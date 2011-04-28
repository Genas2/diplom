# -*- coding: utf-8 -*-

from PyQt4 import QtGui
from PyQt4 import QtCore
from PyQt4.QtCore import SIGNAL, SLOT
from ui_mainwindow import *

import sympy
from equations import Equations

class MainWindow(QtGui.QMainWindow, Ui_MainWindow):

    # Coordinate system mode
    mode = 'cartesian'

    equations = Equations()

    # Cartesian values limits
    min_x = -1
    max_x = 1

    min_y = -1
    max_y = 1

    min_z = -1
    max_z = 1

    # Cartesian spin limits
    lower_x_limit = -100
    upper_x_limit = 100

    lower_y_limit = -100
    upper_y_limit = 100

    lower_z_limit = -100
    upper_z_limit = 100

    # Spherical values limits
    min_r = 0
    max_r = 1

    min_theta = -sympy.pi
    max_theta = sympy.pi

    min_phi = -2 * sympy.pi
    max_phi = 2 * sympy.pi

    # Spherical spin limits
    lower_r_limit = 0
    upper_r_limit = 100

    lower_theta_limit = -2 * sympy.pi
    upper_theta_limit = 2 * sympy.pi

    lower_phi_limit = -2 * sympy.pi
    upper_phi_limit = 2 * sympy.pi

    plot = ''

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

        self.connect(self.spinMin_x_r, SIGNAL('valueChanged(double)'),
                      self, SLOT('changedMinimum_x_r(double)'))
        self.connect(self.spinMin_y_theta, SIGNAL('valueChanged(double)'),
                      self, SLOT('changedMinimum_y_theta(double)'))
        self.connect(self.spinMin_z_phi, SIGNAL('valueChanged(double)'),
                      self, SLOT('changedMinimum_z_phi(double)'))

        self.connect(self.spinMax_x_r, SIGNAL('valueChanged(double)'),
                      self, SLOT('changedMaximum_x_r(double)'))
        self.connect(self.spinMax_y_theta, SIGNAL('valueChanged(double)'),
                      self, SLOT('changedMaximum_y_theta(double)'))
        self.connect(self.spinMax_z_phi, SIGNAL('valueChanged(double)'),
                      self, SLOT('changedMaximum_z_phi(double)'))

        self.connect(self.cmbSysCoord, SIGNAL('currentIndexChanged(const QString)'), 
                      self, SLOT('toggleMode(const QString)'))
        self.connect(self.cmbEquations, SIGNAL('currentIndexChanged(const QString)'), 
                      self, SLOT('setEquation(const QString)'))

        self.connect(self.btnBuildPlot, SIGNAL('clicked()'),
                      self, SLOT('buildPlot()'))
        self.connect(self.btnPreview, SIGNAL('clicked()'),
                      self, SLOT('previewEquation()'))

        mode_id = self.cmbSysCoord.findText(self.mode)
        if mode_id > 0:
            self.cmbSysCoord.setCurrentIndex(mode_id)
            self.cmbSysCoord.emit(SIGNAL('currentIndexChanged(const QString)'), self.mode)

        for line in self.equations.__equations__:
            label, equation, mode = line
            if self.cmbEquations.findText(label):
                self.cmbEquations.addItem(label)

        self.spin_n.setValue(self.equations.n_val)
        self.spin_n.emit(SIGNAL('valueChanged(int)'), self.equations.n_val)
        self.setEquation('Angular part')

    @QtCore.pyqtSlot('int')
    def set_n(self, new_value):
        self.equations.n_val = new_value
        if new_value <= self.equations.l_val:
            self.equations.l_val = new_value - 1
        self.spin_l.setMaximum(new_value - 1)
        self.spin_m.setMaximum(self.equations.l_val)
        self.spin_m.setMinimum(-self.equations.l_val)

    @QtCore.pyqtSlot('int')
    def set_l(self, new_value):
        self.equations.l_val = new_value
        self.spin_m.setMaximum(self.equations.l_val)
        self.spin_m.setMinimum(-self.equations.l_val)

    @QtCore.pyqtSlot('int')
    def set_m(self, new_value):
        self.equations.m_val = new_value

    @QtCore.pyqtSlot('double')
    def changedMinimum_x_r(self, new_value):
        self.spinMax_x_r.setMinimum(new_value)
        if self.mode == 'cartesian':
            self.min_x = new_value
        elif self.mode == 'spherical':
            self.min_r = new_value

    @QtCore.pyqtSlot('double')
    def changedMinimum_y_theta(self, new_value):
        self.spinMax_y_theta.setMinimum(new_value)
        if self.mode == 'cartesian':
            self.min_y = new_value
        elif self.mode == 'spherical':
            self.min_theta = new_value

    @QtCore.pyqtSlot('double')
    def changedMinimum_z_phi(self, new_value):
        self.spinMax_z_phi.setMinimum(new_value)
        if self.mode == 'cartesian':
            self.min_z = new_value
        elif self.mode == 'spherical':
            self.min_phi = new_value

    @QtCore.pyqtSlot('double')
    def changedMaximum_x_r(self, new_value):
        self.spinMin_x_r.setMaximum(new_value)
        if self.mode == 'cartesian':
            self.max_x = new_value
        elif self.mode == 'spherical':
            self.max_r = new_value

    @QtCore.pyqtSlot('double')
    def changedMaximum_y_theta(self, new_value):
        self.spinMin_y_theta.setMaximum(new_value)
        if self.mode == 'cartesian':
            self.max_y = new_value
        elif self.mode == 'spherical':
            self.max_theta = new_value

    @QtCore.pyqtSlot('double')
    def changedMaximum_z_phi(self, new_value):
        self.spinMin_z_phi.setMaximum(new_value)
        if self.mode == 'cartesian':
            self.max_z = new_value
        elif self.mode == 'spherical':
            self.max_phi = new_value

    @QtCore.pyqtSlot()
    def buildPlot(self):
        if type(self.plot).__name__ == 'Plot':
            self.plot.clear()
        else:
            self.plot = sympy.Plot()
        #equation = self.equation(self.equations.l_val, self.equations.m_val)
        if self.mode == 'cartesian':
            self.plot[1] = (self.equation(), [self.equations.x, self.min_x, self.max_x], [self.equations.y, self.min_y, self.max_y])
        elif self.mode == 'spherical':
            if self.equation.__name__ == 'Radial_Part':
                self.plot.axes._label_axes = True
                print(self.max_r)
                self.plot[1] = (self.equation(), [self.equations.r, self.min_r, self.max_r])
            else:
                self.plot[1] = (sympy.trigsimp(self.equation().subs(sympy.sin(self.equations.theta),sympy.sin(self.equations.theta)**2)),
                            [self.equations.phi, self.min_phi, self.max_phi, 35], 
                            [self.equations.theta, self.min_theta, self.max_theta, 35],
                            'mode=spherical')

        self.plot.show()
        return True

    @QtCore.pyqtSlot()
    def previewEquation(self):
        sympy.preview(self.equation())

    @QtCore.pyqtSlot('QString')
    def toggleMode(self, mode):
        if mode == self.mode:
            return 0
        else:
            self.mode = mode

        if mode == 'cartesian':
            # Change coordinates labels
            self.lblMin_x_r.setText('x')
            self.lblMin_y_theta.setText('y')
            self.lblMin_z_phi.setText('z')

            self.lblMax_x_r.setText('x')
            self.lblMax_y_theta.setText('y')
            self.lblMax_z_phi.setText('z')

            # Change spin limits
            self.spinMin_x_r.setMinimum(self.lower_x_limit)
            self.spinMin_y_theta.setMinimum(self.lower_y_limit)
            self.spinMin_z_phi.setMinimum(self.lower_z_limit)

            self.spinMax_x_r.setMaximum(self.upper_x_limit)
            self.spinMax_y_theta.setMaximum(self.upper_y_limit)
            self.spinMax_z_phi.setMaximum(self.upper_z_limit)
            
            # Change values
            self.spinMin_x_r.setValue(self.min_x)
            self.spinMin_y_theta.setValue(self.min_y)
            self.spinMin_z_phi.setValue(self.min_z)

            self.spinMax_x_r.setValue(self.max_x)
            self.spinMax_y_theta.setValue(self.max_y)
            self.spinMax_z_phi.setValue(self.max_z)

        elif mode == 'spherical':
            # Change coordinates labels
            self.lblMin_x_r.setText('r')
            self.lblMin_y_theta.setText(u'θ')
            self.lblMin_z_phi.setText(u'φ')

            self.lblMax_x_r.setText('r')
            self.lblMax_y_theta.setText(u'θ')
            self.lblMax_z_phi.setText(u'φ')

            # Change spin limits
            self.spinMin_x_r.setMinimum(self.lower_r_limit)
            self.spinMin_y_theta.setMinimum(self.lower_theta_limit)
            self.spinMin_z_phi.setMinimum(self.lower_phi_limit)

            self.spinMax_x_r.setMaximum(self.upper_r_limit)
            self.spinMax_y_theta.setMaximum(self.upper_theta_limit)
            self.spinMax_z_phi.setMaximum(self.upper_phi_limit)
            
            # Change values
            self.spinMin_x_r.setValue(self.min_r)
            self.spinMin_y_theta.setValue(self.min_theta)
            self.spinMin_z_phi.setValue(self.min_phi)

            self.spinMax_x_r.setValue(self.max_r)
            self.spinMax_y_theta.setValue(self.max_theta)
            self.spinMax_z_phi.setValue(self.max_phi)
            
    @QtCore.pyqtSlot('QString')
    def setEquation(self, equation):
        for line in self.equations.__equations__:
            label, eq, mode = line
            if label == equation:
                mode_id = self.cmbSysCoord.findText(mode)
                self.cmbSysCoord.setCurrentIndex(mode_id)
                self.equation = eq
            #self.equation = self.equations.Angular_Part

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


