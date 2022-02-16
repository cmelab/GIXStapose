[![build](https://github.com/cmelab/GIXStapose/actions/workflows/build.yml/badge.svg)](https://github.com/cmelab/GIXStapose/actions/workflows/build.yml)
[![pytest](https://github.com/cmelab/GIXStapose/actions/workflows/pytest.yml/badge.svg)](https://github.com/cmelab/GIXStapose/actions/workflows/pytest.yml)
[![codecov](https://codecov.io/gh/cmelab/GIXStapose/branch/master/graph/badge.svg)](https://codecov.io/gh/cmelab/GIXStapose)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/cmelab/GIXStapose/master.svg)](https://results.pre-commit.ci/latest/github/cmelab/GIXStapose/master)


# GIXStapose
GIXStapose is an interactive analysis tool for for studying semi-crystalline and soft materials structures. It enables grazing incidence X-ray scattering (GIXS) patterns to be visualized while interactively rotating chemical structures, especially periodic simulation volumes generated from molecular simulations.

This functionality is useful for interactively identifying  real-space chemical features that correspond to bright diffraction peaks and the rotation matrices that generate them.
As such, this tool has potential to aid in the [reproducible generation of figures that include both GIXS and structural data](http://dx.doi.org/10.1080/08927022.2017.1296958), and it has pedagogical potential for students learning about crystal structures and diffraction.

GIXStapose is made possible by open-source packages, including the high-quality rendering of the [Fresnel ray-tracer](https://fresnel.readthedocs.io/en/stable/), the ability to read in multiple chemical file formats of [MBuild](https://mosdef.org/mbuild/index.html), and numpy's fast Fourier implementations used in interactive [diffraction analysis](https://bitbucket.org/cmelab/cme_utils/src/master/cme_utils/analyze/diffractometer.py), and by funding from the National Science Foundation (#1835593)

![A screen capture of GIXStapose in action](gixstapose/data/screenshot.gif)

### Examples

Click the binder link to see an example using the gixstapose library for generating reproducible figures and scattering patterns:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cmelab/gixstapose/master?urlpath=tree/examples/Figure_Example.ipynb)

### Updates
- Jul 31, 2020: GIXStapose's diffraction code has been added as a module in [Freud](https://freud.readthedocs.io/en/v2.4.1/modules/diffraction.html)!

### Installation
To install GIXStapose you will need the conda package manager (we recommend [Miniconda](https://docs.conda.io/en/latest/miniconda.html))
1. Clone this repo
```
git clone git@github.com:cmelab/GIXStapose.git;
cd GIXStapose
```
2. Using conda, create and activate your environment
```
conda env create -f environment.yml;
conda activate gixstapose
```
3. With the environment active, install this package
```
pip install .
```

### Usage
To run a simple cubic example:
```
gixstapose
```
To load an input file format supported by [MDTraj](http://mdtraj.org/1.8.0/load_functions.html) (e.g., pdb, xml, dcd, xyz, hoomdxml) or a [gsd file](https://gsd.readthedocs.io/en/stable/):
```
gixstapose -i INPUTFILE
```
