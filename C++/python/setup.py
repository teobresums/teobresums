#Always prefer setuptools over distutils
from setuptools import setup, find_packages
from setuptools import Extension
from Cython.Build import cythonize
from setuptools.command.build_ext import build_ext as _build_ext
# To use a consistent encoding
from codecs import open
from os import path
import subprocess
import numpy

# see https://stackoverflow.com/a/21621689/1862861 for why this is here
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

f = open("../config.status","r")

for l in f:
    if 'S["prefix"]' in l:
        library = l.split("=")[1][1:-2]

gsl_include = subprocess.check_output(['gsl-config', '--cflags'])[2:]

ext_modules=[
             Extension("pyTEOBResumS",
                       sources=["pyTEOBResumS/pyTEOBResumS.pyx"],
                       libraries=["m","TEOBResumS","gsl", "gslcblas"], # Unix-like specific
                       include_dirs=[numpy.get_include(),"../src",gsl_include],
                       language='c++',
                       extra_compile_args=["-Ofast"],
                       extra_link_args=["-L"+library+"/lib/"]
                       )
             ]

setup(
      name = "pyTEOBResumS",
      version = '1.0',
      description = 'pyTEOBResumS: Spinning Tidal EOB model',
      author = 'Walter Del Pozzo',
      author_email='walter.delpozzo@unipi.it',
      ext_modules = cythonize(ext_modules),
      include_dirs=[numpy.get_include(),"../src",gsl_include],
      license='MIT',
      cmdclass = {'build_ext': build_ext},
      classifiers =[
                    # How mature is this project? Common values are
                    #   3 - Alpha
                    #   4 - Beta
                    #   5 - Production/Stable
                    'Development Status :: 4 - Beta',
                    
                    # Indicate who your project is intended for
                    #'Intended Audience :: Developers',
                    #'Topic :: Data Analysis :: Bayesian Inference',
                    
                    # Pick your license as you wish (should match "license" above)
                    'License :: OSI Approved :: MIT License',
                    
                    # Specify the Python versions you support here. In particular, ensure
                    # that you indicate whether you support Python 2, Python 3 or both.
                    'Programming Language :: Python :: 2',
                    'Programming Language :: Python :: 2.6',
                    'Programming Language :: Python :: 2.7',
                    'Programming Language :: Python :: 3',
                    'Programming Language :: Python :: 3.2',
                    'Programming Language :: Python :: 3.3',
                    'Programming Language :: Python :: 3.4',
                    'Programming Language :: Python :: 3.5',
                    ],
      
      keywords='gravitational waves',
      
      packages=['pyTEOBResumS'],
      install_requires=['numpy'],
      setup_requires=['numpy'],
      package_data={"": ['*.pyx', '*.pxd']},
      # To provide executable scripts, use entry points in preference to the
      # "scripts" keyword. Entry points provide cross-platform support and allow
      # pip to create the appropriate form of executable for the target platform.
      entry_points={
      #    'console_scripts':['sample=sample:main',
      #        ],
      }
      )
