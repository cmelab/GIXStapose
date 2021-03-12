from setuptools import setup

VERSION = None

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    with open(os.path.join(here, NAME, "__version__.py")) as f:
        exec(f.read(), about)
else:
    about["__version__"] = VERSION

setup(name='gixstapose',
      version=about["__version__"],
      description='An interactive structure viewer alongside its diffraction pattern',
      url='https://github.com/cmelab/GIXStapose',
      author='Jenny Fothergill',
      author_email='jennyfothergill@boisestate.edu',
      license='BSD-3',
      packages=['gixstapose'],
      package_dir={'gixstapose': 'gixstapose'},
      package_data={'gixstapose': ['data']},
      zip_safe=False,
      entry_points = {
          'console_scripts': ['gixstapose=gixstapose.main:main'],
          }
      )
