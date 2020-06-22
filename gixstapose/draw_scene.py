import random
import warnings

import matplotlib.cm
import numpy as np
import gsd
import gsd.pygsd
import gsd.hoomd

import fresnel
import freud
import mbuild as mb
import PIL

cpk_colors = {
    "H": fresnel.color.linear([1.00, 1.00, 1.00]),  # white
    "C": fresnel.color.linear([0.30, 0.30, 0.30]),  # grey
    "N": fresnel.color.linear([0.13, 0.20, 1.00]),  # dark blue
    "O": fresnel.color.linear([1.00, 0.13, 0.00]),  # red
    "F": fresnel.color.linear([0.12, 0.94, 0.12]),  # green
    "Cl": fresnel.color.linear([0.12, 0.94, 0.12]),  # green
    "Br": fresnel.color.linear([0.60, 0.13, 0.00]),  # dark red
    "I": fresnel.color.linear([0.40, 0.00, 0.73]),  # dark violet
    "He": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "Ne": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "Ar": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "Xe": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "Kr": fresnel.color.linear([0.00, 1.00, 1.00]),  # cyan
    "P": fresnel.color.linear([1.00, 0.60, 0.00]),  # orange
    "S": fresnel.color.linear([1.00, 0.90, 0.13]),  # yellow
    "B": fresnel.color.linear([1.00, 0.67, 0.47]),  # peach
    "Li": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Na": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "K": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Rb": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Cs": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Fr": fresnel.color.linear([0.47, 0.00, 1.00]),  # violet
    "Be": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Mg": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Ca": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Sr": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Ba": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Ra": fresnel.color.linear([0.00, 0.47, 0.00]),  # dark green
    "Ti": fresnel.color.linear([0.60, 0.60, 0.60]),  # grey
    "Fe": fresnel.color.linear([0.87, 0.47, 0.00]),  # dark orange
    "default": fresnel.color.linear([0.87, 0.47, 1.00]),  # pink
}

bsu_colors = [
        fresnel.color.linear([0.00, 0.20, 0.63]), # blue
        fresnel.color.linear([0.84, 0.30, 0.04]), # orange
        fresnel.color.linear([0.00, 0.23, 0.44]), # blue
        fresnel.color.linear([0.00, 0.42, 0.65]), # blue
        fresnel.color.linear([0.00, 0.45, 0.81]), # blue
        fresnel.color.linear([0.25, 0.38, 0.60]), # blue
        fresnel.color.linear([0.17, 0.99, 1.00]), # blue
        fresnel.color.linear([1.00, 0.39, 0.03]), # orange
        fresnel.color.linear([1.00, 0.40, 0.00]), # orange
        ]

# Made space to add more later
radii_dict = {"H": 0.05, "default": 0.06}


def from_gsd(gsdfile, frame=-1, coords_only=False, scale=1.0):
    """
    Given a trajectory gsd file creates an mbuild.Compound.
    If there are multiple separate molecules, they are returned
    as one compound.

    Parameters
    ----------
    gsdfile : str, filename
    frame : int, frame number (default -1)
    coords_only : bool (default False)
        If True, return compound with no bonds
    scale : float, scaling factor multiplied to coordinates (default 1.0)

    Returns
    -------
    mbuild.Compound
    """
    f = gsd.pygsd.GSDFile(open(gsdfile, 'rb'))
    t = gsd.hoomd.HOOMDTrajectory(f)

    snap = t[frame]
    bond_array = snap.bonds.group
    n_atoms = snap.particles.N

    # Add particles
    comp = mb.Compound()
    comp.box = snap.configuration.box[:3] * scale
    for i in range(n_atoms):
        name = snap.particles.types[snap.particles.typeid[i]]
        xyz = snap.particles.position[i] * scale
        charge = snap.particles.charge[i]

        atom = mb.Particle(name=name, pos=xyz, charge=charge)
        comp.add(atom, label=str(i))

    if not coords_only:
        # Add bonds
        for i in range(bond_array.shape[0]):
            atom1 = int(bond_array[i][0])
            atom2 = int(bond_array[i][1])
            comp.add_bond([comp[atom1], comp[atom2]])
    return comp

def distance(pos1, pos2):
    """
    Calculates euclidean distance between two points.

    Parameters
    ----------
    pos1, pos2 : ((3,) numpy.ndarray), xyz coordinates
        (2D also works)

    Returns
    -------
    float distance
    """
    return np.linalg.norm(pos1 - pos2)


def mb_to_freud_box(box):
    """
    Convert an mbuild box object to a freud box object
    These sites are helpful as reference:
    http://gisaxs.com/index.php/Unit_cell
    https://hoomd-blue.readthedocs.io/en/stable/box.html

    Parameters
    ----------
    box : mbuild.box.Box()

    Returns
    -------
    freud.box.Box()
    """
    Lx = box.lengths[0]
    Ly = box.lengths[1]
    Lz = box.lengths[2]
    alpha = box.angles[0]
    beta = box.angles[1]
    gamma = box.angles[2]

    frac = (
        np.cos(np.radians(alpha)) - np.cos(np.radians(beta)) * np.cos(np.radians(gamma))
    ) / np.sin(np.radians(gamma))
    c = np.sqrt(1 - np.cos(np.radians(beta)) ** 2 - frac ** 2)

    xy = np.cos(np.radians(gamma)) / np.sin(np.radians(gamma))
    xz = frac / c
    yz = np.cos(np.radians(beta)) / c

    box_list = list(box.maxs) + [xy, yz, xz]
    return freud.box.Box(*box_list)


def visualize(comp, color="cpk", scale=1.0, box=None):
    """
    Visualize an mbuild Compound using fresnel.

    Parameters
    ----------
    comp : (mbuild.Compound), compound to visualize
    color : ("cpk", "bsu", or the name of a matplotlib colormap), color scheme to use (default "cpk")
    scale : (float), scaling factor for the particle, bond, and box radii (default 1.0)
    box : (mb.Box), box object for the structure. If no box is provided,
        comp.boundingbox is used

    Returns
    -------
    fresnel.Scene
    """

    ## Extract some info about the compound
    N = comp.n_particles
    particle_names = [p.name for p in comp.particles()]

    N_bonds = comp.n_bonds
    if N_bonds > 0:
        # all_bonds.shape is (nbond, 2 ends, xyz)
        all_bonds = np.stack([np.stack((i[0].pos, i[1].pos)) for i in comp.bonds()])

    color_array = np.empty((N, 3), dtype="float64")
    if color == "cpk":
        # Populate the color_array with colors based on particle name
        # -- if name is not defined in the dictionary, use pink (the default)
        for i, n in enumerate(particle_names):
            try:
                color_array[i, :] = cpk_colors[n]
            except KeyError:
                color_array[i, :] = cpk_colors["default"]
    elif color == "bsu":
        # Populate the color array with the brand standard bsu colors
        # https://www.boisestate.edu/communicationsandmarketing/brand-standards/colors/
        # if there are more unique particle names than colors, colors will be reused
        unique_names = list(set(particle_names))
        for i, n in enumerate(particle_names):
            color_array[i, :] = bsu_colors[uniq_atoms.index(n) % len(bsu_colors)]
    else:
        # Populate the color_array with colors based on particle name
        # choose colors evenly distributed through a matplotlib colormap
        try:
            cmap = matplotlib.cm.get_cmap(name=color)
        except ValueError:
            print(
                "The 'color' argument takes either 'cpk', 'bsu', or the name of a matplotlib colormap."
            )
            raise
        mapper = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(vmin=0, vmax=1, clip=True), cmap=cmap
        )
        particle_types = list(set(particle_names))
        N_types = len(particle_types)
        v = np.linspace(0, 1, N_types)
        # Color by typeid
        type_ids = np.array([particle_types.index(i) for i in particle_names])
        for i in range(N_types):
            color_array[type_ids == i] = fresnel.color.linear(mapper.to_rgba(v)[i])

    # Make an array of the radii based on particle name
    # -- if name is not defined in the dictionary, use default
    rad_array = np.empty((N), dtype="float64")
    for i, n in enumerate(particle_names):
        try:
            rad_array[i] = radii_dict[n] * scale
        except KeyError:
            rad_array[i] = radii_dict["default"] * scale

    ## Start building the fresnel scene
    scene = fresnel.Scene()

    # Spheres for every particle in the system
    geometry = fresnel.geometry.Sphere(scene, N=N)
    geometry.position[:] = comp.xyz
    geometry.material = fresnel.material.Material(roughness=1.0)
    geometry.outline_width = 0.01 * scale

    # use color instead of material.color
    geometry.material.primitive_color_mix = 1.0
    geometry.color[:] = color_array

    # resize radii
    geometry.radius[:] = rad_array

    # bonds
    if N_bonds > 0:
        bonds = fresnel.geometry.Cylinder(scene, N=N_bonds)
        bonds.material = fresnel.material.Material(roughness=0.5)
        bonds.outline_width = 0.01 * scale

        # bonds are white
        bond_colors = np.ones((N_bonds, 3), dtype="float64")

        bonds.material.primitive_color_mix = 1.0
        bonds.points[:] = all_bonds

        bonds.color[:] = np.stack(
            [fresnel.color.linear(bond_colors), fresnel.color.linear(bond_colors)], axis=1
        )
        bonds.radius[:] = [0.03 * scale] * N_bonds

    if box:
        # Use comp.box, unless it does not exist, then use comp.boundingbox
        try:
            freud_box = mb_to_freud_box(box)
        except AttributeError:
            freud_box = mb_to_freud_box(comp.boundingbox)
        # Create box in fresnel
        fresnel.geometry.Box(scene, freud_box, box_radius=0.008 * scale)
    return scene


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


class CG(mb.Compound):
    def __init__(self):
        super(CG, self).__init__()
        np.random.seed(42)
        for i in range(ord("A"), ord("Z") + 1):
            random.seed(i)
            N = random.randint(1, 4)
            for n in range(N):
                self.add(
                    mb.Particle(name=f"_{chr(i)}", pos=np.random.random(3) * 2 - 1)
                )

        for i in range(self.n_particles - 1):
            for j in range(i + 1, self.n_particles):
                dist = distance(self.xyz[i], self.xyz[j])
                if dist < 0.3:
                    self.add_bond((self[i], self[j]))
