import fresnel
import mbuild as mb
from PySide2 import QtCore, QtGui, QtWidgets

import interact
from draw_scene import Methane, visualize

# Build example scene
methane = Methane()
methane.box = mb.Box(lengths=[0.5, 0.5, 0.5])

scene = visualize(methane, show_box=True)

class _Camera(QtWidgets.QLabel):
    def __init__(self, text="", *args, **kwargs):
        super(_Camera, self).__init__(*args, **kwargs)
        self.setText(text)

    def _trigger_refresh(self):
        self.update()


class Diffract(QtWidgets.QWidget):
    def __init__(self, scene=scene, *args, **kwargs):
        super().__init__(*args, **kwargs)

        layout = QtWidgets.QHBoxLayout()
        self._view = interact.SceneView(scene)

        self._camera = _Camera(str(self._view.scene.camera))

        self._view.changeEvent.connect(self._camera._trigger_refresh)
        layout.addWidget(self._view)

        layout.addWidget(self._camera)
        self.setLayout(layout)
