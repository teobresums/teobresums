"""
Setup script for CPython_test.c
Compile with 
$ python CPython_test_setup.py build_ext --inplace
"""

from distutils.core import setup, Extension
import numpy
setup(ext_modules=[
    Extension("foo_module", ["CPython_test.c"],
              include_dirs=[numpy.get_include()],
              extra_compile_args=["-std=c99"]),],)



