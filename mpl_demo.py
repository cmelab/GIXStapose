import sys

import numpy as np

from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
if is_pyqt5():
    print("using pyqt5")
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
else:
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)

        # Create the QDial widget and set up defaults.
        # - we provide accessors on this class to override.
        self._dial = QtWidgets.QDial()
        self._dial.setNotchesVisible(True)
        self._dial.setWrapping(False)
        self._dial.setMinimum(0)
        self._dial.setMaximum(10)
        layout.addWidget(self._dial)

        dynamic_canvas = FigureCanvas(Figure(figsize=(5, 3)))
        layout.addWidget(dynamic_canvas)
        #self.addToolBar(QtCore.Qt.BottomToolBarArea,
        #                NavigationToolbar(dynamic_canvas, self))

        self._dynamic_ax = dynamic_canvas.figure.subplots()
        t = np.linspace(0, 10, 101)
        self._dynamic_ax.plot(t, np.sin(t))
        self._dial.valueChanged.connect(self._update_canvas)

    def _update_canvas(self):
        self._dynamic_ax.clear()
        val = self._dial.value()
        # Shift the sinusoid as a function of time.
        t = np.linspace(0, 10, 101)
        self._dynamic_ax.plot(t, np.sin(t))
        self._dynamic_ax.plot(val, np.sin(val),"ro")
        self._dynamic_ax.figure.canvas.draw()


if __name__ == "__main__":
    qapp = QtWidgets.QApplication(sys.argv)
    app = ApplicationWindow()
    app.show()
    qapp.exec_()

