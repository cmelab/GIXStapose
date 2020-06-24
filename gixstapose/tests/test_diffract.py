from pathlib import Path
from tempfile import NamedTemporaryFile
from PIL import Image

import matplotlib.pyplot as plt
import numpy as np

from gixstapose.draw_scene import compound_load
from gixstapose.diffractometer import Diffractometer, camera_to_rot
from gixstapose.main import camera_from_pos

path = str(Path(__file__).parent.parent.resolve())
temp = NamedTemporaryFile(suffix=".png")

def test_diffract():
    d = Diffractometer()
    inputfile = path + "/data/sc10.pdb"
    compound = compound_load(inputfile)
    d.load(compound.xyz,compound.boundingbox.lengths)
    rot = camera_to_rot(
            camera_from_pos((1,0,0))
            )
    dp = d.diffract(rot.T)
    plt.imsave(temp.name, dp, cmap="jet")
    dpim = Image.open(temp.name)
    dparr = np.asarray(dpim)
    im = Image.open(path + "/data/sc10_camera100.png")
    imarr = np.asarray(im)
    assert np.allclose(dparr,imarr)


