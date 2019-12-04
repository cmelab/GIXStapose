import numpy as np

from matplotlib.backends.qt_compat import is_pyqt5 # QtCore, QtWidgets, is_pyqt5
from PySide2 import QtCore, QtWidgets
if is_pyqt5():
    print("using pyqt5")
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
else:
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure


import fresnel
import mbuild as mb

import interact
from draw_scene import Methane, visualize

# Build example scene
methane = Methane()
methane.box = mb.Box(lengths=[0.5, 0.5, 0.5])

scene = visualize(methane, show_box=True)

class ApplicationWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        layout = QtWidgets.QHBoxLayout()

        # Add the SceneView widget
        self._view = interact.SceneView(scene)
        layout.addWidget(self._view)

        #self._view.c.update_camera.connect(self._update_camera)

        #dynamic_canvas = FigureCanvas(Figure(figsize=(5, 3)))
        #layout.addWidget(dynamic_canvas)
        ##self.addToolBar(QtCore.Qt.BottomToolBarArea,
        ##                NavigationToolbar(dynamic_canvas, self))

        #self._dynamic_ax = dynamic_canvas.figure.subplots()
        #t = np.linspace(0, 10, 101)
        #self._dynamic_ax.plot(t, np.sin(t))
        #self._dial.valueChanged.connect(self._update_canvas)

        self._label = QtWidgets.QLabel()
        self._label.setText(str(self._view.scene.camera))
        layout.addWidget(self._label)
        #self._dial.valueChanged.connect(self._update_text)

        self.setLayout(layout)

    #def _update_canvas(self):
    #    self._dynamic_ax.clear()
    #    val = self._dial.value()
    #    # Shift the sinusoid as a function of time.
    #    t = np.linspace(0, 10, 101)
    #    self._dynamic_ax.plot(t, np.sin(t))
    #    self._dynamic_ax.plot(val, np.sin(val),"ro")
    #    self._dynamic_ax.figure.canvas.draw()

    def _update_camera(self):
        self._label.clear()
        val = self._view.scene.camera
        # display the camera value
        self._label.setText(str(val))


if __name__ == "__main__":
    qapp = QtWidgets.QApplication([])

    app = ApplicationWindow()
    app.show()

    qapp.exec_()

