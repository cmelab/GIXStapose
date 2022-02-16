import fresnel
import pytest
from base_test import BaseTest

from gixstapose.draw_scene import get_scene


class TestScene(BaseTest):
    def test_scene_gsd(self, methanegsd):
        scene, info = get_scene(methanegsd)
        assert type(scene) is type(fresnel.Scene())

        scene, info = get_scene(methanegsd, show_bonds=True)
        assert type(scene) is type(fresnel.Scene())

    def test_scene_comp(self, methanemol2):
        scene, info = get_scene(methanemol2)
        assert type(scene) is type(fresnel.Scene())

        scene, info = get_scene(methanemol2, show_bonds=True)
        assert type(scene) is type(fresnel.Scene())

    def test_scene_colors(self, methanegsd):
        scene, info = get_scene(methanegsd, color="bsu")
        assert type(scene) is type(fresnel.Scene())

        scene, info = get_scene(methanegsd, color={"C": "grey"})
        assert type(scene) is type(fresnel.Scene())

        scene, info = get_scene(methanegsd, color="jet")
        assert type(scene) is type(fresnel.Scene())

        with pytest.raises(ValueError):
            get_scene(methanegsd, color="heck")
