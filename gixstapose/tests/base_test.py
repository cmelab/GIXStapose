import os
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from gixstapose.diffractometer import camera_to_rot
from gixstapose.draw_scene import get_info
from gixstapose.main import camera_from_pos

data_dir = str(Path(__file__).parent.parent.resolve()) + "/data/"


class BaseTest:
    @pytest.fixture
    def sc10(self):
        return os.path.join(data_dir, "sc10.pdb")

    @pytest.fixture
    def positions_and_box(self, sc10):
        info = get_info(sc10)
        return info["positions"], info["box"][:3]

    @pytest.fixture
    def camera100(self):
        return camera_from_pos((1, 0, 0))

    @pytest.fixture
    def rot100(self, camera100):
        return camera_to_rot(camera100).T

    @pytest.fixture
    def imarray100(self):
        im = Image.open(os.path.join(data_dir, "sc10_camera100.png"))
        return np.asarray(im)

    @pytest.fixture
    def methanegsd(self):
        return os.path.join(data_dir, "methane.gsd")

    @pytest.fixture
    def methanemol2(self):
        return os.path.join(data_dir, "methane.mol2")
