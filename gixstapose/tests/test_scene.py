import os

import fresnel
import pytest

from gixstapose.draw_scene import get_scene


test_dir = os.path.dirname(__file__)
data_dir = os.path.join(test_dir, "../data")


def test_scene_gsd():
    methanegsd = os.path.join(data_dir, "methane.gsd")
    scene, info = get_scene(methanegsd)
    assert type(scene) is type(fresnel.Scene())

    scene, info = get_scene(methanegsd, show_bonds=True)
    assert type(scene) is type(fresnel.Scene())

def test_scene_comp():
    methanemol2 = os.path.join(data_dir, "methane.mol2")
    scene, info = get_scene(methanemol2)
    assert type(scene) is type(fresnel.Scene())

    scene, info = get_scene(methanemol2, show_bonds=True)
    assert type(scene) is type(fresnel.Scene())

def test_scene_colors():
    methanegsd = os.path.join(data_dir, "methane.gsd")
    scene, info = get_scene(methanegsd, color="bsu")
    assert type(scene) is type(fresnel.Scene())

    scene, info = get_scene(methanegsd, color={"C": "grey"})
    assert type(scene) is type(fresnel.Scene())

    scene, info = get_scene(methanegsd, color="jet")
    assert type(scene) is type(fresnel.Scene())

    with pytest.raises(ValueError):
        get_scene(methanegsd, color="heck")
