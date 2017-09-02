from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import subprocess

f = open("../config.status","r")

for l in f:
    if 'S["prefix"]' in l:
        library = l.split("=")[1][1:-2]

gsl_include = subprocess.check_output(['gsl-config', '--cflags'])[2:]

ext_modules=[
             Extension("pyTEOBResum",
                       sources=["pyTEOBResum.pyx"],
                       libraries=["m","TEOBResum","gsl", "gslcblas"], # Unix-like specific
                       include_dirs=[numpy.get_include(),"../src",gsl_include],
                       language='c++',
                       extra_compile_args=["-O3"],
                       extra_link_args=["-L"+library+"/lib/"]
                       )
             ]

setup(
      name = "pyTEOBResum",
      ext_modules = cythonize(ext_modules),
      include_dirs=[numpy.get_include(),"../src",gsl_include]
      )
