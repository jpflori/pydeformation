#!/usr/bin/env python
## -*- encoding: utf-8 -*-

import os
import sys
from setuptools import setup
from codecs import open # To open the README file with proper encoding
from setuptools.command.test import test as TestCommand # for tests
from setuptools.extension import Extension
from Cython.Build import cythonize
from sage.env import sage_include_directories, SAGE_LOCAL

# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename,  encoding='utf-8') as f:
        return f.read()

# For the tests
class SageTest(TestCommand):
    def run_tests(self):
        errno = os.system("sage -t --force-lib pydeformation")
        if errno != 0:
            sys.exit(1)

if not os.path.isfile(os.path.join(SAGE_LOCAL, "include", "deformation", "deformation.h")):
    print("The deformation library is not installed.")
    sys.exit(1)

cythonize_dir = "build"

kwds = {"include_dirs": sage_include_directories()}

extensions = [
    Extension("pydeformation.deformation", ["pydeformation/deformation.pyx"], **kwds)
]

setup(
    name="pydeformation",
    author="Jean-Pierre Flori",
    author_email="sage-devel@googlegroups.com",
    url="https://github.com/jpflori/pydeformation",
    license="GNU General Public License, version 3 or later",
    description="Wrapper for deformation library by Sebastian Pancratz",
    long_description = readfile("README.rst"), # get the long description from the README
    version = readfile("VERSION"), # the VERSION file is shared with the documentation
    classifiers=[
      # How mature is this project? Common values are
      #   3 - Alpha
      #   4 - Beta
      #   5 - Production/Stable
      'Development Status :: 4 - Beta',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Mathematics',
      'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
      'Programming Language :: Python :: 2.7',
    ], # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords = "sagemath deformation",
    setup_requires=["cython", "sagemath"],
    install_requires=["sagemath"],
    ext_modules = cythonize(extensions),
    packages=["pydeformation"],
    include_package_data = True,
    cmdclass = {'test': SageTest} # adding a special setup command for tests
)
