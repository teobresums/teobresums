"""
Setup script for TEOBResumSWrap.c based on distutils
"""

from setuptools import setup, Extension
import numpy
import glob
import os
import re

# We use env variables to set macros (debug, verbose, user_funcs, etc)
DEBUG_MODE = os.environ.get("TEOB_DEBUG", "0") == "1"
VERBOSE_MODE = os.environ.get("TEOB_VERBOSE", "0") == "1"
USE_TIDAL = os.environ.get("TEOB_USEBTIDALPOTENTIAL", "1") == "1"
USER_FUNCS = os.environ.get("TEOB_USER_FUNCTIONS", "0") == "1"

# Read version from pyproject.toml
ini = open("pyproject.toml").read()
vrs = r"^version = ['\"]([^'\"]*)['\"]"
mo = re.search(vrs, ini, re.M)
version = mo.group(1)

# Our module
eob_sources = ["TEOBResumSWrap.c"]

# Sources
if not os.path.exists("lib"):
    os.symlink("../C/src", "lib")
if not os.path.exists("PyREADME.md"):
    os.symlink("../README.md", "PyREADME.md")

eob_sources += glob.glob("lib/*.c")

# Include dirs
incdirs = [numpy.get_include(), "lib"]

# Set macros and compiler args based on env variables
macros = [("HAVE_GSL", 1)]
extra_compile_args = ["-w", "-std=c99", "-ffast-math"]

if DEBUG_MODE:
    print("Compiling in DEBUG mode")
    macros.append(("DEBUG", 1))
    extra_compile_args.append("-g")
else:
    macros.append(("DEBUG", 0))
    extra_compile_args.append("-O3")

if VERBOSE_MODE:
    print("Compiling in VERBOSE mode")
    macros.append(("VERBOSE", 1))
else:
    macros.append(("VERBOSE", 0))

# Handle Tidal state
if USE_TIDAL:
    macros.append(("USEBTIDALPOTENTIAL", 1))
else:
    macros.append(("USEBTIDALPOTENTIAL", 0))

# Handle user functions state
if USER_FUNCS:
    macros.append(("USERFUNCS", 1))
    user_src = "../C/user/TEOBResumSUser.c"
    # also include user funcs. This is only for local
    # source compilation
    if os.path.exists(user_src):
        eob_sources.append(user_src)
        incdirs.append("../C/user")

# Libraries
libs = ["gsl", "gslcblas", "m"]

# copy src files to this repo

pyprofit_ext = Extension(
    "EOBRun_module",
    depends=glob.glob("lib/*.h"),
    language="c",
    define_macros=macros,
    sources=eob_sources,
    include_dirs=incdirs,
    libraries=libs,
    extra_compile_args=extra_compile_args,
    # Binds calls between functions defined in this module directly, instead
    # of routing every cross-file call (eob_metric_s, eob_dyn_s_GS, ...) through
    # the PLT to allow for symbol interposition that never happens for a
    # single self-contained extension module. ~19% faster in-process on the
    # M=2 Dali BBH benchmark, bit-identical output (see Utils/Profiling/README.md).
    extra_link_args=["-Wl,-Bsymbolic-functions"],
)

setup(
    name="teobresums",
    version=version,
    description="TEOBResumS: Effective-One-Body model for generic compact binary coalescences",
    author="The TEOBResumS team",
    url="https://bitbucket.org/eob_ihes/teobresums/src/master/",
    ext_modules=[pyprofit_ext],
    install_requires=["setuptools>42", "numpy>=1.18.0"],
)
