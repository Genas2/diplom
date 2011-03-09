#!/usr/bin/python2

import sys, errno, os
from PyQt4 import QtGui, uic

ui_file = 'mainwindow.ui'

# Compiled to python ui file
uic_file = 'ui_mainwindow.py'
uic_fd = open(uic_file, 'w')

try:
    uic.compileUi(open(ui_file, 'r'), uic_fd)
except IOError:
    print(sys.argv[0] + ": " + ui_file + ' ' + os.strerror(errno.ENOENT))
    exit(errno.ENOENT)

uic_fd.close()

from mainwindow import MainWindow

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
