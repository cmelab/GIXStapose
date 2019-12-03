import fresnel
import fresnel.interact
import mbuild as mb

from draw_scene import Methane, visualize

# build Scene
methane = Methane()
methane.box = mb.Box(lengths=[0.5,0.5,0.5])

scene = visualize(methane, show_box=True)


view = fresnel.interact.SceneView(scene)
view.show()
fresnel.interact.app.exec_()
