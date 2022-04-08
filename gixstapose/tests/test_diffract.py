import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pytest
from base_test import BaseTest
from PIL import Image

from gixstapose.diffractometer import Diffractometer


class Test_Diffractometer(BaseTest):
    def test_diffract(self, positions_and_box, rot100, imarray100, tmp_path):
        d_file = tmp_path / "dp.png"
        d = Diffractometer()
        d.load(*positions_and_box)
        dp = d.diffract(rot100)
        plt.imsave(d_file, dp, cmap="jet")
        dpim = Image.open(d_file)
        dparr = np.asarray(dpim)
        assert np.allclose(dparr, imarray100)

    def test_diffract_plot_camera(self, positions_and_box, camera100):
        d = Diffractometer(length_scale=1.0)
        d.load(*positions_and_box)
        d.diffract_from_camera(camera100)
        fig, ax = d.plot()
        assert isinstance(ax, plt.Axes)
        assert isinstance(fig, plt.Figure)
        assert all(np.isclose(ax.get_xlim(), (-405, 405), atol=1))

    def test_diffract_plot_rot(self, positions_and_box, rot100):
        d = Diffractometer(length_scale=2.0)
        d.load(*positions_and_box)
        d.diffract(rot100)
        fig, ax = d.plot()
        assert isinstance(ax, plt.Axes)
        assert isinstance(fig, plt.Figure)
        assert all(np.isclose(ax.get_xlim(), (-202, 202), atol=1))

    def test_diffract_plot_raises(self, positions_and_box):
        with pytest.raises(ValueError):
            d = Diffractometer()
            d.plot()
        with pytest.raises(ValueError):
            d = Diffractometer()
            d.load(*positions_and_box)
            d.plot()
