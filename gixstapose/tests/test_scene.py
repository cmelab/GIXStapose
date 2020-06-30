import fresnel
import mbuild as mb

from gixstapose.draw_scene import create_scene

class Methane(mb.Compound):
    def __init__(self):
        super(Methane, self).__init__()
        carbon = mb.Particle(name="C")
        self.add(carbon, label="C[$]")

        hydrogen = mb.Particle(name="H", pos=[0.1, 0, -0.07])
        self.add(hydrogen, label="HC[$]")

        self.add_bond((self[0], self["HC"][0]))

        self.add(mb.Particle(name="H", pos=[-0.1, 0, -0.07]), label="HC[$]")
        self.add(mb.Particle(name="H", pos=[0, 0.1, 0.07]), label="HC[$]")
        self.add(mb.Particle(name="H", pos=[0, -0.1, 0.07]), label="HC[$]")

        self.add_bond((self[0], self["HC"][1]))
        self.add_bond((self[0], self["HC"][2]))
        self.add_bond((self[0], self["HC"][3]))

def test_visualise():
    methane = Methane()
    scene = create_scene(methane)
    assert type(scene) is type(fresnel.Scene())

