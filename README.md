# diffract
An interactive structure viewer with alongside its simulated diffraction pattern
This code uses [fresnel](https://fresnel.readthedocs.io/en/stable/), [mbuild](https://mosdef.org/mbuild/index.html), and [cme_utils/diffractometer](https://bitbucket.org/cmelab/cme_utils/src/master/cme_utils/analyze/diffractometer.py)

To install:
1. Create environment and install necessary packages
```
conda create -n diffract python=3.7
conda activate diffract
conda install -c conda-forge -c omnia -c mosdef pillow numpy matplotlib mbuild fresnel pyside2 freud py3dmol openbabel jupyterlab;
```
2. Clone and install the cme_utils package
```
git clone git@bitbucket.org:cmelab/cme_utils.git
cd cme_utils
pip install -e .
```
