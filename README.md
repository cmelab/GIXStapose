# diffract
An interactive structure viewer alongside its simulated diffraction pattern
This code uses [fresnel](https://fresnel.readthedocs.io/en/stable/), [mbuild](https://mosdef.org/mbuild/index.html), and [cme_utils/diffractometer](https://bitbucket.org/cmelab/cme_utils/src/master/cme_utils/analyze/diffractometer.py)

![A screen capture of diffract in action](screenshot.gif)

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

To run:
```
python main.py
```
to run a simple-cubic example or
```
python main.py -i INPUTFILE
```
to load an input file format supported by [MDTraj](http://mdtraj.org/1.8.0/load_functions.html) (e.g., pdb, xml, dcd, xyz, hoomdxml) or a [gsd file](https://gsd.readthedocs.io/en/stable/).

