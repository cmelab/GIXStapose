import numpy as np

from matplotlib.backends.qt_compat import is_pyqt5 # QtCore, QtWidgets, is_pyqt5
from PySide2 import QtCore, QtWidgets
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
else:
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure


import fresnel
import mbuild as mb

from diffractometer import Diffractometer
import interact
from draw_scene import visualize #Methane

# Build example scene
pdbname = "sc10"
dirname = "gixs_data"

box = mb.Box(np.array([10,10,10])/10)
pdb = mb.load(f"{dirname}/{pdbname}.pdb")
pdb.box = box

scene = visualize(pdb, show_box=True)

class ApplicationWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.title = "diffractometer"
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)

        self.createGridLayout()

        windowlayout = QtWidgets.QVBoxLayout()
        windowlayout.addWidget(self.horizontalGroupBox)

        self._label = QtWidgets.QLabel()
        self._label.setText(self.getCameraText(self._view.scene.camera))
        windowlayout.addWidget(self._label)

        self.setLayout(windowlayout)

        self.show()


        #dynamic_canvas = FigureCanvas(Figure(figsize=(5, 3)))
        #layout.addWidget(dynamic_canvas)
        ##self.addToolBar(QtCore.Qt.BottomToolBarArea,
        ##                NavigationToolbar(dynamic_canvas, self))

        #self._dynamic_ax = dynamic_canvas.figure.subplots()
        #t = np.linspace(0, 10, 101)
        #self._dynamic_ax.plot(t, np.sin(t))
        #self._dial.valueChanged.connect(self._update_canvas)

        #self._dial.valueChanged.connect(self._update_text)

    def getCameraText(self, camera):
        pos = camera.position
        look = camera.look_at
        text =  "".join([
            "Camera\n",
            "   position : {0:.3f} {1:.3f} {2:.3f}\n".format(pos[0], pos[1], pos[2]),
            "   look at :  {0:.3f} {1:.3f} {2:.3f}\n".format(look[0], look[1], look[2]),
            "   up :       {0:.3f} {1:.3f} {2:.3f}\n".format(
                camera.up[0], camera.up[1], camera.up[2]
                ),
            "   height :   {0:.3f}".format(camera.height)
            ])
        return text

    def createGridLayout(self):
        self.horizontalGroupBox = QtWidgets.QGroupBox()
        layout = QtWidgets.QGridLayout()
        #layout.setColumnStretch(1, 4)
        #layout.setColumnStretch(2, 4)

        # Add the SceneView widget
        self._view = interact.SceneView(scene)
        self._view.c.update_camera.connect(self._update_camera)
        layout.addWidget(self._view,0,0)

        # Add the diffraction widget
        self._diffract = QtWidgets.QLabel()
        self._diffract.setText("diffraction pattern")
        layout.addWidget(self._diffract,0,1)

        self.horizontalGroupBox.setLayout(layout)

    #def _update_canvas(self):
    #    self._dynamic_ax.clear()
    #    val = self._dial.value()
    #    # Shift the sinusoid as a function of time.
    #    t = np.linspace(0, 10, 101)
    #    self._dynamic_ax.plot(t, np.sin(t))
    #    self._dynamic_ax.plot(val, np.sin(val),"ro")
    #    self._dynamic_ax.figure.canvas.draw()

    def _update_camera(self, value):
        self._label.clear()
        # display the camera value
        self._label.setText(self.getCameraText(value))

    def diffract(self, rot=None):
        pass



if __name__ == "__main__":
    qapp = QtWidgets.QApplication([])

    app = ApplicationWindow()
    app.show()

    qapp.exec_()

