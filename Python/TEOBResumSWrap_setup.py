"""
Setup script for TEOBResumSWrap.c based on distutils
Compile with 
$ python TEOBResumSWrap_setup.py build_ext --inplace
"""

from distutils.core import setup, Extension
import distutils.ccompiler
import numpy
import glob


# Our module
eob_sources = ['TEOBResumSWrap.c']

# Sources
eob_sources += glob.glob('../C/src/*.c')

# Include dirs
incdirs = [numpy.get_include(), '/usr/include/']

# Libraries
libs = ['gsl', 'gslcblas', 'm', 'config']


pyprofit_ext = Extension('EOBRun_module',
                         depends=glob.glob('src/*.h'),
                         language='c99',
                         define_macros = [('HAVE_GSL',1)],
                         sources = eob_sources,
                         include_dirs = incdirs,
                         libraries = libs,
                         extra_compile_args=["-w", "-std=c99"] #to ignore all warnings (do not use when debugging)
                       )

setup(ext_modules=[pyprofit_ext],)
