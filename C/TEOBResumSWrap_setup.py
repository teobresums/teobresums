"""
Setup script for TEOBResumSWrap.c
Compile with 
$ python TEOBResumSWrap_setup.py build_ext --inplace

Inspired by https://github.com/ICRAR/pyprofit/issues/4
"""

from distutils.core import setup, Extension
import distutils.ccompiler
import numpy
import glob
#setup(ext_modules=[Extension("EOBRunTD_module", ["TEOBResumSWrap.c", "TEOBResumSDynamics.c", "TEOBResumSEOBRun.c", "TEOBResumSFits.c", "TEOBResumSFlux.c", "TEOBResumSInitialConditions.c", "TEOBResumSMetric.c", "TEOBResumSOMP.c", "TEOBResumSPars.c", "TEOBResumSPostAdiab.c", "TEOBResumSUtils.c", "TEOBResumSWaveform.c", "TEOBResumS.c"], include_dirs=[numpy.get_include()]),],)

# Our module
eob_sources = ['TEOBResumSWrap.c']

# sources
eob_sources += glob.glob('src/*.c')

# include dirs
incdirs = [numpy.get_include(), '/usr/include/']

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



