import os
import random

import fresnel
import gsd.hoomd
import matplotlib.cm
from matplotlib import colors as mplcolors
import mbuild as mb
import numpy as np
import PIL

from gixstapose.color_dicts import cpk_colors, bsu_colors, radii_dict


def get_scene(
    inputfile, frame=-1, color="cpk", scale=1.0, show_bonds=False, scene=None
):
    """Loads an input file into a fresnel.Scene.

    Parameters
    ----------
    inputfile: str,
        Path to input file
    frame: int,
        Which frame to load if inputfile is a trajectory. Supports negative
        indexing (default -1)
    color : str, default "cpk"
        Color scheme to use
        ("cpk", "bsu", name of a matplotlib colormap, or a custom dictionary)
    scale : float, default 1.0
        Scaling factor for the particle radii, bond and box lengths
    show_bonds : bool, default False
        Whether to show bonds
    scene : fresnel.Scene, default None
        Existing scene to add to. If None is provided, a new scene is created.

    Returns
    -------
    fresnel.Scene
    """
    info = get_info(inputfile, frame, show_bonds)

    return create_scene(info, color, scale, show_bonds, scene), info


def create_scene(
    info, color="cpk", scale=1.0, show_bonds=False, init_scene=None
):
    """Create a fresnel.Scene object.

    Adds geometries for particles, bonds, and box (or boundingbox).

    Parameters
    ----------
    info : list
        List containing N, types, typeids, positions, N_bonds, bonds, box
    color : str, default "cpk"
        Color scheme to use
        ("cpk", "bsu", name of a matplotlib colormap, or a custom dictionary)
    scale : float, default 1.0
        Scaling factor for the particle radii, bond and box lengths
    show_bonds : bool, default False
        Whether to show bonds
    init_scene : fresnel.Scene, default None
        Existing scene to add to. If None is provided, a new scene is created.

    Returns
    -------
    fresnel.Scene
    """

    N, types, typeids, positions, N_bonds, bonds, box = info

    color_array = np.empty((N, 3), dtype="float64")
    if type(color) is dict:
        # Populate the color_array with colors based on particle name
        # if name is not defined in the dictionary, try using the cpk dictionary
        for i, n in enumerate(typeids):
            try:
                ncolor = color[types[n]]
                color_array[i, :] = fresnel.color.linear(
                    mplcolors.to_rgba(ncolor)
                )
            except KeyError:
                try:
                    color_array[i, :] = cpk_colors[types[n]]
                except KeyError:
                    color_array[i, :] = cpk_colors["default"]
    elif color == "cpk":
        # Populate the color_array with colors based on particle name
        # -- if name is not defined in the dictionary, use pink (the default)
        for i, n in enumerate(typeids):
            try:
                color_array[i, :] = cpk_colors[types[n]]
            except KeyError:
                color_array[i, :] = cpk_colors["default"]
    elif color == "bsu":
        # Populate the color array with the brand standard bsu colors
        # https://www.boisestate.edu/communicationsandmarketing/
        # brand-standards/colors/
        # if there are more unique particle names than colors,
        # colors will be reused
        for i, n in enumerate(typeids):
            color_array[i, :] = bsu_colors[n % len(bsu_colors)]
    else:
        # Populate the color_array with colors based on particle name
        # choose colors evenly distributed through a matplotlib colormap
        try:
            cmap = matplotlib.cm.get_cmap(name=color)
        except ValueError:
            print(
                "The 'color' argument takes either 'cpk', 'bsu', or the name "
                "of a matplotlib colormap."
            )
            raise
        mapper = matplotlib.cm.ScalarMappable(
            norm=mplcolors.Normalize(vmin=0, vmax=1, clip=True), cmap=cmap
        )
        N_types = len(types)
        v = np.linspace(0, 1, N_types)
        # Color by typeid
        for i in range(N_types):
            color_array[typeids == i] = fresnel.color.linear(
                mapper.to_rgba(v)[i]
            )

    # Make an array of the radii based on particle name
    # -- if name is not defined in the dictionary, use default
    rad_array = np.empty((N), dtype="float64")
    for i, n in enumerate(typeids):
        try:
            rad_array[i] = radii_dict[types[n]] * scale
        except KeyError:
            rad_array[i] = radii_dict["default"] * scale

    ## Start building the fresnel scene
    if init_scene is None:
        scene = fresnel.Scene()
    else:
        scene = init_scene

    # Spheres for every particle in the system
    geometry = fresnel.geometry.Sphere(scene, N=N)
    geometry.position[:] = positions
    geometry.material = fresnel.material.Material(roughness=1.0)
    geometry.outline_width = 0.01 * scale

    # use color instead of material.color
    geometry.material.primitive_color_mix = 1.0
    geometry.color[:] = color_array

    # resize radii
    geometry.radius[:] = rad_array

    # bonds
    if N_bonds > 0 and show_bonds:
        bond_cyls = fresnel.geometry.Cylinder(scene, N=N_bonds)
        bond_cyls.material = fresnel.material.Material(roughness=0.5)
        bond_cyls.outline_width = 0.01 * scale

        # bonds are white
        bond_colors = np.ones((N_bonds, 3), dtype="float64")

        bond_cyls.material.primitive_color_mix = 1.0
        bond_cyls.points[:] = bonds

        bond_cyls.color[:] = np.stack(
            [
                fresnel.color.linear(bond_colors),
                fresnel.color.linear(bond_colors)
            ], axis=1
        )
        bond_cyls.radius[:] = [0.03 * scale] * N_bonds

    if init_scene is None:
        # Create box in fresnel
        fresnel.geometry.Box(scene, box, box_radius=0.008 * scale)

        # Set the initial camera position
        max_dist = np.max(positions) - np.min(positions)
        scene.camera.height = 1.5 * max_dist
        scene.camera.position = [max_dist, max_dist, max_dist]
    return scene


def get_comp_info(comp, show_bonds):
    N = comp.n_particles
    names = [p.name for p in comp.particles()]
    types = list(set(names))
    typeids = np.array([types.index(i) for i in names])
    positions = comp.xyz

    N_bonds = comp.n_bonds
    if N_bonds > 0 and show_bonds:
        # bonds.shape is (nbond, 2 ends, xyz)
        bonds = np.stack(
            [np.stack((i[0].pos, i[1].pos)) for i in comp.bonds()]
        )
    else:
        bonds = None

    # Use comp.box, unless it does not exist, then use comp.boundingbox
    if comp.box is not None:
        box = comp.box
    else:
        box = comp.get_boundingbox()
    box = np.hstack([box.lengths, box.tilt_factors])

    return [N, types, typeids, positions, N_bonds, bonds, box]


def get_gsd_info(gsdfile, frame, show_bonds):
    with gsd.hoomd.open(gsdfile) as f:
        snap = f[frame]

    N = snap.particles.N
    types = snap.particles.types
    typeids = snap.particles.typeid
    positions = snap.particles.position

    N_bonds = snap.bonds.N
    if N_bonds > 0 and show_bonds:
        # bonds.shape is (nbond, 2 ends, xyz)
        bonds = np.stack(
            [(positions[i],positions[j]) for (i,j) in snap.bonds.group]
        )
    else:
        bonds = None

    box = snap.configuration.box
    return [N, types, typeids, positions, N_bonds, bonds, box]


def get_info(inputfile, frame=-1, show_bonds=False):
    name, extension = os.path.splitext(inputfile)
    if extension == ".gsd":
        # info contains N, types, typeids, positions, N_bonds, bonds, box
        info = get_gsd_info(inputfile, frame, show_bonds)
    else:
        comp = mb.load(inputfile)
        # info contains N, types, typeids, positions, N_bonds, bonds, box
        info = get_comp_info(comp, show_bonds)
    return info
