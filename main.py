import numpy as np

from PySide2 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import fresnel
import mbuild as mb

from diffractometer import Diffractometer, camera_to_rot
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

        # Initialize the diffractometer
        self.d = Diffractometer()
        self.d.load(pdb.xyz, box.maxs)

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

        # Add the SceneView widget
        self._view = interact.SceneView(scene)
        self._view.c.update_camera.connect(self._update_camera)
        layout.addWidget(self._view,0,0)

        # Add the diffraction widget
        dynamic_canvas = FigureCanvas(Figure(figsize=(15,15)))
        layout.addWidget(dynamic_canvas,0,1)
        self._diffract_ax = dynamic_canvas.figure.add_subplot(111)
        self._diffract_ax.axis("off")
        self._diffract_ax.set_axis_off()
        self.plot_diffract(self._view.scene.camera)

        self.horizontalGroupBox.setLayout(layout)

    def _update_camera(self, value):
        self._label.clear()
        # display the camera value
        self._label.setText(self.getCameraText(value))
        self.plot_diffract(value)

    def plot_diffract(self, camera):
        self._diffract_ax.clear()
        self._diffract_ax.axis("off")
        rot = camera_to_rot(camera)

        # diffraction pattern
        dp = self.d.diffract(rot.T)
        self._diffract_ax.imshow(dp, cmap="jet")
        self._diffract_ax.figure.canvas.draw()





if __name__ == "__main__":
    qapp = QtWidgets.QApplication([])

    app = ApplicationWindow()
    app.show()

    qapp.exec_()

