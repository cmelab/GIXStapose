"""GIXStapose QT application."""

import argparse
import os
from pathlib import Path

import fresnel
import matplotlib.pyplot as plt
import numpy as np
import PIL
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# PySide2 must be imported before matplotlib -- isort will switch these
from PySide2.QtCore import QSize, Qt
from PySide2.QtWidgets import (
    QAction,
    QApplication,
    QGridLayout,
    QGroupBox,
    QLabel,
    QMainWindow,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from gixstapose.diffractometer import Diffractometer
from gixstapose.draw_scene import get_scene


class ApplicationWindow(QMainWindow):  # pragma: no cover
    """Main class to hold all GIXStapose application."""

    def __init__(self, inputfile, frame):
        super().__init__()

        self.title = "GIXStapose"

        if inputfile is None:
            print("No input provided, showing simple cubic example")
            path = Path(__file__).parent / "data/sc10.pdb"
            inputfile = str(path.resolve())
        self.basename = os.path.basename(inputfile).split(".")[0]
        self.render_counter = 0
        self.diffract_counter = 0

        self.init_diffractometer(inputfile, frame)

        self.initUI()

    def init_diffractometer(self, inputfile, frame):
        """Initialize the diffractometer."""
        self.scene, info = get_scene(inputfile, frame)
        mode = self.scene.device.mode
        self.title = f"GIXStapose ({mode} mode)"

        self.d = Diffractometer()
        self.d.load(info["positions"], info["box"][:3])

    def initUI(self):
        """Initialize the user interface."""
        self.setWindowTitle(self.title)

        # Menubar
        menubar = self.menuBar()
        filemenu = menubar.addMenu("File")
        render = QAction("Render Scene", self)
        export = QAction("Export Diffraction Pattern", self)
        print_camera = QAction("Print Camera", self)
        filemenu.addAction(render)
        filemenu.addAction(export)
        filemenu.addAction(print_camera)
        filemenu.triggered[QAction].connect(self.processtrigger)

        self.main = QWidget()
        self.setCentralWidget(self.main)

        # creates 'top' and 'bot' horizontalgroupbox grid layout objects
        self.createGridLayout()

        windowlayout = QVBoxLayout()
        windowlayout.addWidget(self.tophorizontalGroupBox)
        windowlayout.addWidget(self.bothorizontalGroupBox)
        self.main.setLayout(windowlayout)

        self.show()

    def createGridLayout(self):
        """Create the grid layout to position widgets."""
        from fresnel import interact

        # Top grid with sceneview and diffraction pattern
        self.tophorizontalGroupBox = QGroupBox()
        toplayout = QGridLayout()

        # Add the SceneView widget
        self.view = interact.SceneView(self.scene)
        self.view.rendering.connect(self.update_camera)
        toplayout.addWidget(self.view, 0, 0)

        # Add the diffraction widget
        dynamic_canvas = FigureCanvas(Figure(figsize=(15, 15)))
        toplayout.addWidget(dynamic_canvas, 0, 1)
        self.diffract_ax = dynamic_canvas.figure.add_subplot(111)
        self.diffract_ax.axis("off")
        self.diffract_ax.set_axis_off()
        self.plot_diffract(self.view.scene.camera)

        self.tophorizontalGroupBox.setLayout(toplayout)

        # Bottom grid with camera, buttons, zoom, sigma
        self.bothorizontalGroupBox = QGroupBox()
        botlayout = QGridLayout()

        # Camera printout
        self.label = QLabel()
        self.label.setText(self.camera_text(self.view.scene.camera))
        # widget, row, col, rowspan, colspan
        botlayout.addWidget(self.label, 0, 0, 2, 1)

        # Buttons
        self.button100 = QPushButton("100")
        self.button100.setMaximumSize(QSize(100, 40))
        botlayout.addWidget(self.button100, 0, 1, 2, 1)

        self.button110 = QPushButton("110")
        self.button110.setMaximumSize(QSize(100, 40))
        botlayout.addWidget(self.button110, 0, 2, 2, 1)

        self.button111 = QPushButton("111")
        self.button111.setMaximumSize(QSize(100, 40))
        botlayout.addWidget(self.button111, 0, 3, 2, 1)

        # Connect buttons to moving the camera
        # thanks to this wonderful answer
        # https://stackoverflow.com/a/57167056/11969403
        self.button100.clicked.connect(lambda: self.move_camera((1, 0, 0)))
        self.button110.clicked.connect(lambda: self.move_camera((1, 1, 0)))
        self.button111.clicked.connect(lambda: self.move_camera((1, 1, 1)))

        # Add space between buttons and slider
        botlayout.setColumnMinimumWidth(4, 350)

        # Zoom slider
        zlabel = QLabel("Zoom")
        zlabel.setAlignment(Qt.AlignCenter)
        botlayout.addWidget(zlabel, 0, 5)

        self.zooms = [i for i in range(1, self.d.N) if self.d.N % i == 0]
        self.zoomslider = QSlider(Qt.Horizontal)
        self.zoomslider.setMinimum(0)
        self.zoomslider.setMaximum(len(self.zooms) - 1)
        self.zoomslider.setValue(self.zooms.index(self.d.zoom))
        self.zoomslider.valueChanged.connect(self.change_zoom)
        self.zoomslider.setMaximumSize(QSize(600, 30))
        botlayout.addWidget(self.zoomslider, 1, 5)

        botlayout.setColumnMinimumWidth(6, 50)

        self.bothorizontalGroupBox.setLayout(botlayout)

    def processtrigger(self, q):
        """Process the file menu triggers."""
        if q.text() == "Render Scene":
            print("Rendering...")
            output = fresnel.pathtrace(
                self.view.scene, light_samples=40, w=600, h=600
            )
            filename = f"{self.basename}_scene{self.render_counter}.png"

            image = PIL.Image.fromarray(output[:], mode="RGBA")
            image.save(filename, dpi=(300, 300))

            old_camera = self.view.scene.camera
            print(f"Rendered {filename}")
            self.render_counter += 1

        elif q.text() == "Export Diffraction Pattern":
            print("Saving diffraction pattern...")
            filename = f"{self.basename}_dp{self.diffract_counter}.png"

            plt.imsave(filename, self.dp, cmap="jet")
            print(f"Diffraction pattern saved as {filename}")
            self.diffract_counter += 1

        elif q.text() == "Print Camera":
            print("Current camera is:\n")
            print(self.camera_text(self.view.scene.camera))

    def change_zoom(self):
        """Change the zoom of the diffractometer."""
        self.d.zoom = self.zooms[self.zoomslider.value()]
        self.plot_diffract(self.view.scene.camera)

    def camera_text(self, camera):
        """Convert a fresnel.Camera object to a readable string."""
        pos = camera.position
        look = camera.look_at
        up = camera.up
        text = "".join(
            [
                "Camera\n",
                f" position = [{pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f}],\n",
                f" look_at =  [{look[0]:.3f}, {look[1]:.3f}, {look[2]:.3f}],\n",
                f" up =       [{up[0]:.3f}, {up[1]:.3f}, {up[2]:.3f}],\n",
                f" height =   {camera.height:.3f}",
            ]
        )
        return text

    def update_camera(self, camera):
        """Update the camera."""
        self.label.clear()
        # display the camera value
        self.label.setText(self.camera_text(camera))
        self.plot_diffract(camera)

    def plot_diffract(self, camera):
        """Plot the difraction pattern."""
        self.diffract_ax.clear()
        self.diffract_ax.axis("off")

        # diffraction pattern
        self.dp = self.d.diffract_from_camera(camera)
        self.diffract_ax.imshow(self.dp, cmap="jet")
        self.diffract_ax.figure.canvas.draw()
        self.repaint()

    def move_camera(self, pos):
        """Move the camera to a position."""
        self.view.scene.camera = camera_from_pos(pos)
        # self.repaint()
        self.view._start_rendering()
        self.view.update()


def camera_from_pos(pos):
    """Create a new camera instance at position."""
    camera = fresnel.camera.Orthographic(
        position=pos, look_at=(0, 0, 0), up=(0, 0, 1), height=1.5
    )

    return camera


def main():  # pragma: no cover
    """Initialize application window from command line."""
    parser = argparse.ArgumentParser(
        description="Provide a chemical input file"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="an input file, accepted types: mol2, pdb, xyz, gsd",
    )
    parser.add_argument(
        "-t",
        "--frame",
        type=int,
        default="-1",
        help="if trajectory file is given, which frame to diffract",
    )
    args = parser.parse_args()

    qapp = QApplication.instance()
    if qapp == None:
        qapp = QApplication([])

    app = ApplicationWindow(args.input, args.frame)
    app.show()

    qapp.exec_()


if __name__ == "__main__":
    main()
