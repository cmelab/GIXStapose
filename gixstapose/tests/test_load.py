from pathlib import Path

import mbuild as mb

from gixstapose.draw_scene import compound_load

path = str(Path(__file__).parent.parent.resolve())

def test_load_pdb():
    inputfile = path + "/data/sc10.pdb"
    print(inputfile)
    compound = compound_load(inputfile)
    assert type(compound) is type(mb.Compound())

def test_load_gsd():
    inputfile = path + "/data/sc10.gsd"
    print(inputfile)
    compound = compound_load(inputfile)
    assert type(compound) is type(mb.Compound())
