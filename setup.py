from setuptools import setup

setup(name='gixstapose',
      version='0.0',
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
