#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
opj = os.path.join

from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
from sage.env import sage_include_directories

cythonize_dir = "build"

kwds = {"include_dirs": sage_include_directories()}

extensions = [
    Extension("pydeformation.deformation", ["src/pydeformation/deformation.pyx"], **kwds)
]

setup(
    name="pydeformation",
    author=u"Edgar Costa, Jean-Pierre Flori",
    author_email="sage-devel@googlegroups.com",
    url="https://github.com/jpflori/pydeformation",
    license="GNU Lesser General Public License, version 3 or later",
    description="Wrapper for deformation library",
    setup_requires=["cython"],
    ext_modules = cythonize(extensions),
    packages=["pydeformation"],
    package_dir={"pydeformation": opj("src", "pydeformation")},
    package_data={"pydeformation": ["*.pxd"]},
)
