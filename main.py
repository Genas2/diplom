#!/usr/bin/python2

import sys, errno, os
from PyQt4 import QtGui, uic
from mainwindow import MainWindow

ui_file = 'mainwindow.ui'

try:
    uic.compileUi(open(ui_file), open('ui_mainwindow.py', 'w'))
except IOError:
    print(sys.argv[0] + ": " + ui_file + ' ' + os.strerror(errno.ENOENT))
    exit(errno.ENOENT)

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
