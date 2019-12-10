import argparse

import fresnel
import mbuild as mb
import numpy as np
# PySide2 must be imported before matplotlib -- isort will switch these
from PySide2 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import interact
from diffractometer import Diffractometer, camera_to_rot
from draw_scene import visualize, from_gsd



class ApplicationWindow(QtWidgets.QWidget):
    def __init__(self, inputfile, frame):
        super().__init__()
        self.title = "diffractometer"

        # Initialize the diffractometer
        self.init_diffractometer(inputfile, frame)

        self.initUI()

    def init_diffractometer(self, inputfile, frame):
        if inputfile is None:
            print("no input provided, showing simple cubic example")
            inputfile = f"example_inputs/sc10.pdb"
        try:
            compound = from_gsd(inputfile, frame=frame)
        except RuntimeError:
            compound = mb.load(inputfile)

        self.scene = visualize(compound)

        self.d = Diffractometer()
        box = compound.boundingbox.maxs - compound.boundingbox.mins
        self.d.load(compound.xyz, box)

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
        text = "".join(
            [
                "Camera\n",
                "   position : {0:.3f} {1:.3f} {2:.3f}\n".format(
                    pos[0], pos[1], pos[2]
                ),
                "   look at :  {0:.3f} {1:.3f} {2:.3f}\n".format(
                    look[0], look[1], look[2]
                ),
                "   up :       {0:.3f} {1:.3f} {2:.3f}\n".format(
                    camera.up[0], camera.up[1], camera.up[2]
                ),
                "   height :   {0:.3f}".format(camera.height),
            ]
        )
        return text

    def createGridLayout(self):
        self.horizontalGroupBox = QtWidgets.QGroupBox()
        layout = QtWidgets.QGridLayout()

        # Add the SceneView widget
        self._view = interact.SceneView(self.scene)
        self._view.c.update_camera.connect(self._update_camera)
        layout.addWidget(self._view, 0, 0)

        # Add the diffraction widget
        dynamic_canvas = FigureCanvas(Figure(figsize=(15, 15)))
        layout.addWidget(dynamic_canvas, 0, 1)
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
    parser = argparse.ArgumentParser(description="Provide a chemical input file")
    parser.add_argument("-i", "--input", type=str,
        help="an input file, accepted types: mol2, pdb, xyz, gsd")
    parser.add_argument("-t", "--frame", type=int, default="-1",
        help="if trajectory file is given, which frame to diffract (default -1 or last frame)")
    args = parser.parse_args()

    qapp = QtWidgets.QApplication([])

    app = ApplicationWindow(args.input, args.frame)
    app.show()

    qapp.exec_()
