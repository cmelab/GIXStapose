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
        return data_dir + "sc10.pdb"

    @pytest.fixture
    def positions_and_box(self, sc10):
        _, _, _, positions, _, _, box = get_info(sc10)
        return positions, box[:3]

    @pytest.fixture
    def camera100(self):
        return camera_from_pos((1, 0, 0))

    @pytest.fixture
    def rot100(self, camera100):
        return camera_to_rot(camera100).T

    @pytest.fixture
    def imarray100(self):
        im = Image.open(data_dir + "sc10_camera100.png")
        return np.asarray(im)
