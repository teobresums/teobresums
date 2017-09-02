from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

f = open("../config.status","r")

for l in f:
    if 'S["prefix"]' in l:
        library = l.split("=")[1][1:-2]


ext_modules=[
             Extension("pyTEOBResum",
                       sources=["pyTEOBResum.pyx"],
                       libraries=["m","TEOBResum"], # Unix-like specific
                       include_dirs=[numpy.get_include(),"../src"],
                       language='c++',
                       extra_compile_args=["-O3"],
                       extra_link_args=["-L"+library+"/lib/"]
                       )
             ]

setup(
      name = "pyTEOBResum",
      ext_modules = cythonize(ext_modules),
      include_dirs=[numpy.get_include(),"../src"]
      )
