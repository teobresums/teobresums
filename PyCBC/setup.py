#!/usr/bin/env python

"""
setup.py file for TEOBResumS pycbc waveform plugin package
"""

from setuptools import Extension, setup, Command
from setuptools import find_packages

NAME = 'teobresums-pycbc'
VERSION = 'v2'

setup (
    name = NAME,
    version = VERSION,
    description = 'PyCBC pluging for the native, stand-alone and actively developed implementation of the TEOBResumS approximant',
    long_description = open('descr.rst').read(),
    author = 'The TEOBResumS team',
    author_email = 'sebastiano.bernuzzi@uni-jena.de',
    url = 'https://bitbucket.org/eob_ihes/teobresums/wiki/Home',
    download_url = 'https://bitbucket.org/eob_ihes/teobresums/src/master/',
    keywords = ['pycbc', 'signal processing', 'gravitational waves'],
    install_requires = ['pycbc'],
    py_modules = ['teobresums'],
    entry_points = {"pycbc.waveform.td":"teobresums = teobresums:teobresums_td",
                    "pycbc.waveform.fd":"teobresums_fd = teobresums:teobresums_fd"},
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
