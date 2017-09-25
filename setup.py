#!/usr/bin/env python
req=['nose','numpy','matplotlib','h5py','pandas','scipy']
# %%
try:
    import conda.cli
    conda.cli.main('install',*req)
except Exception as e:
    import pip
    pip.main(['install'] + req)
# %%
from setuptools import setup

setup(name='gnsseclipse',
      description='utilities for the ionospheric remote sensing of Eclipse 2017 via GNSS',
      author='Sebastijan Mrak',
      url='https://github.com/aldebaran1/gnsseclipse',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 3 - Alpha',
      'Topic :: Scientific/Engineering :: Atmospheric Science',
      'Programming Language :: Python',
      'Programming Language :: Python :: 3',
      ],
      install_requires=req,
      packages=['gnsseclipse']
)