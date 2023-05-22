"""
Setup script for TEOBResumSWrap.c based on distutils
Compile with 
$ python setup.py build_ext --inplace
"""

from distutils.core import setup, Extension
import numpy; import glob; import os; import re

# Read version from pyproject.toml
ini = open('pyproject.toml').read()
vrs = r"^version = ['\"]([^'\"]*)['\"]"
mo  = re.search(vrs, ini, re.M)
version = mo.group(1)

# Our module
eob_sources = ['TEOBResumSWrap.c']

# Sources
if not os.path.isdir('lib')         : os.symlink('../C/src', 'lib')
if not os.path.exists("PyREADME.md"): os.symlink("../README.md", "PyREADME.md")

eob_sources += glob.glob('lib/*.c')

# Include dirs
incdirs = [numpy.get_include(), '/usr/include/', 'lib']

# Libraries
libs = ['gsl', 'gslcblas', 'm']

# copy src files to this repo

pyprofit_ext = Extension('EOBRun_module',
                         depends=glob.glob('lib/*.h'),
                         language='c99',
                         define_macros = [('HAVE_GSL',1)],
                         sources = eob_sources,
                         include_dirs = incdirs,
                         libraries = libs,
                         extra_compile_args=["-w", "-std=c99", "-DDEBUG=0", "-DVERBOSE=0"]#, "-fcommon"
                       )

setup(  name='teobresums',
        version=version,
        description='TEOBResumS: Effective-One-Body model for generic compact binary coalescences',
        author = 'The TEOBResumS team',
        url = 'https://bitbucket.org/eob_ihes/teobresums/src/master/',
        ext_modules=[pyprofit_ext],
        install_requires=['setuptools>42', 'numpy>=1.18.0'],
    )
