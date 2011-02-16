from PyQt4 import QtCore, QtGui
from ui_mainwindow import *
#from world import WorldArea

class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self, parent)

        self.setupUi(self)
        self.saveAsActs = []

        #self.world = WorldArea()
        #self.setCentralWidget(self.world)

        self.createActions()
        self.createMenus()

        self.setWindowTitle(self.tr("Life"))
        self.resize(500, 500)

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


