from PySide2 import QtCore, QtGui, QtWidgets

from diffract import Diffract


app = QtWidgets.QApplication([])

diffract = Diffract()
diffract.show()

app.exec_()
