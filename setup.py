#!/usr/bin/env python3
from sys import platform
import glob
import pkg_resources
import os
import numpy as np
from setuptools import setup, find_packages, Extension
from sysconfig import get_paths
from Cython.Build import cythonize
from Cython.Compiler import Options

Options.language_level = 3

__INCLUDE = get_paths()['include'].split(os.sep)[:-1]
__INCLUDE = (os.sep).join(__INCLUDE)

if platform == "win32":
    __CANTERA_OBJ = pkg_resources.resource_filename('cantera', '_cantera.*lib')
else:
    __CANTERA_OBJ = pkg_resources.resource_filename('cantera', '_cantera.*so')
__CANTERA_OBJ = glob.glob(__CANTERA_OBJ)[0]

__CANTERA_DEP = pkg_resources.resource_filename('cantera', 'interrupts.py')


def readme():
    with open('README.md') as f:
        return f.read()


exts = [Extension("ctmicro._ctmicro",
                  ["ctmicro/_ctmicro.pyx",
                   "ctmicro/Profile.cpp",
                   "ctmicro/ChannelFlow.cpp", ],
                  extra_objects=[__CANTERA_OBJ],
                  include_dirs=[np.get_include(), __INCLUDE],
                  depends=[__CANTERA_DEP],), ]

setup(name="ctmicro",
      description='Non-adiabatic combustion based on Cantera',
      long_description=readme(),
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: BSD 3',
          'Programming Language :: Python :: 3',
      ],
      packages=find_packages(),
      install_requires=['cantera>=2.5.0a4',
                        'cython>=0.29.7',
                        'pandas'],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      ext_modules=cythonize(exts))
