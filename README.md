# diffract
interactive structure viewer with diffraction

To install:

```
conda create -n diffract python=3.7
conda activate diffract
conda install -c conda-forge -c omnia -c mosdef pillow pytest sphinx sphinx_rtd_theme nbsphinx numpy matplotlib black isort jupyterlab mbuild fresnel pyside2 freud py3dmol openbabel;
jupyter labextension install @ryantam626/jupyterlab_code_formatter;
pip install jupyterlab_code_formatter;
jupyter serverextension enable --py jupyterlab_code_formatter;
```

clone the cme_utils package

```
git clone git@bitbucket.org:cmelab/cme_utils.git
cd cme_utils
pip install -e .
```
