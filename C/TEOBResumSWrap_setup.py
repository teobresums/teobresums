"""
Setup script for CPython_test.c
Compile with 
$ python CPython_test_setup.py build_ext --inplace
"""

from distutils.core import setup, Extension
import distutils.ccompiler
import numpy
import glob
#setup(ext_modules=[Extension("EOBRunTD_module", ["TEOBResumSWrap.c", "TEOBResumSDynamics.c", "TEOBResumSEOBRun.c", "TEOBResumSFits.c", "TEOBResumSFlux.c", "TEOBResumSInitialConditions.c", "TEOBResumSMetric.c", "TEOBResumSOMP.c", "TEOBResumSPars.c", "TEOBResumSPostAdiab.c", "TEOBResumSUtils.c", "TEOBResumSWaveform.c", "TEOBResumS.c"], include_dirs=[numpy.get_include()]),],)

#def has_gsl():
#    compiler = distutils.ccompiler.new_compiler()
#    compiler.add_include_dir('/usr/include') #RAPH
#    compiler.add_library_dir('/usr/local/lib') #RAPH
#    return compiler.has_function('gsl_sf_gamma', libraries=['gsl', 'gslcblas'])
#    return compiler.has_function('gsl_sf_gamma', libraries=['/usr/local/lib/gsl', '/usr/local/lib/gslcblas'])

# Our module
eob_sources = ['TEOBResumSWrap.c']

# sources
eob_sources += glob.glob('src/*.c')

# include dirs
#incdirs = ['libprofit', '/usr/local/include'] #RAPH
#incdirs = ['libprofit', '/usr/local/include']
incdirs = [numpy.get_include(), '/usr/include/']

# gsl libs
#if not has_gsl():
#    print("\n\nNo GSL installation found on your system. Install the GSL development package and try again\n\n")
#    sys.exit(1)

#libs = ['gsl', 'gslcblas']
libs = ['gsl', 'gslcblas', 'm', 'config']

pyprofit_ext = Extension('EOBRunTD_module',
                       depends=glob.glob('src/*.h'),
                       language='c99',
                       define_macros = [('HAVE_GSL',1)],
                       sources = eob_sources,
                       include_dirs = incdirs,
                       libraries = libs
                       #extra_compile_args=["-lgsl -lgslcblas -lm -lconfig"]
                       )

setup(ext_modules=[pyprofit_ext],)



