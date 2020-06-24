![CI](https://github.com/cmelab/GIXStapose/workflows/CI/badge.svg)
# GIXStapose
GIXStapose is a new interactive analysis tool for for studying semi-crystalline and soft materials structures. It enables grazing incidence X-ray scattering (GIXS) patterns to be visualized while interactively rotating chemical structures, especially periodic simulation volumes generated from molecular simulations.

This functionality is useful for interactively identifying  real-space chemical features that correspond to bright diffraction peaks and the rotation matrices that generate them.
As such, this tool has potential to aid in the [reproducible generation of figures that include both GIXS and structural data](http://dx.doi.org/10.1080/08927022.2017.1296958), and it has pedagogical potential for students learning about crystal structures and diffraction.

GIXStapose is made possible by open-source packages, including the high-quality rendering of the [Fresnel ray-tracer](https://fresnel.readthedocs.io/en/stable/), the ability to read in multiple chemical file formats of [MBuild](https://mosdef.org/mbuild/index.html), and numpy's fast Fourier implementations used in [interactive diffraction analysis](https://bitbucket.org/cmelab/cme_utils/src/master/cme_utils/analyze/diffractometer.py), and by funding from the National Science Foundation (#1835593)

![A screen capture of GIXStapose in action](gixstapose/data/screenshot.gif)

To install GIXStapose you will need the conda package manager (we recommend [Miniconda](https://docs.conda.io/en/latest/miniconda.html))
1. Using conda, create and activate your environment
```
conda env create -f environment.yml;
conda activate gixstapose
```
2. Clone this repo
```
git clone git@github.com:cmelab/GIXStapose.git;
cd GIXStapose
```
3. Install this package with pip
```
pip install -e .
```

To run a simple cubic example:
```
gixstapose
```
to load an input file format supported by [MDTraj](http://mdtraj.org/1.8.0/load_functions.html) (e.g., pdb, xml, dcd, xyz, hoomdxml) or a [gsd file](https://gsd.readthedocs.io/en/stable/):
```
gixstapose -i INPUTFILE
```

Related citations:

[Jones, M. L., & Jankowski, E. (2017). Computationally connecting organic photovoltaic performance to atomistic arrangements and bulk morphology. Molecular Simulation, 43(10–11), 1–18. https://doi.org/10.1080/08927022.2017.1296958](http://dx.doi.org/10.1080/08927022.2017.1296958)

[Miller, E. D., Jones, M. L., Henry, M. M., Chery, P., Miller, K., & Jankowski, E. (2018). Optimization and Validation of Efficient Models for Predicting Polythiophene Self-Assembly. Polymers, 10(12), 1305. https://doi.org/10.3390/polym10121305](https://doi.org/10.3390/polym10121305)
