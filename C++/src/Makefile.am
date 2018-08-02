AM_LDFLAGS=`gsl-config --libs` -lgsl -lgslcblas -lm
ACLOCAL_AMFLAGS='-I m4 `gsl-config --cflags`'

lib_LTLIBRARIES = libTEOBResumS.la
libTEOBResumS_la_SOURCES =           \
    TEOBResumSFlux.cpp               \
    TEOBResumSMetric.cpp             \
    TEOBResumSUtils.cpp              \
    TEOBResumS.cpp                   \
    TEOBResumSHlm.cpp                \
    TEOBResumSWaveform.cpp           \
    TEOBResumSInitialConditions.cpp  \
    TEOBResumSRingdown.cpp

bin_PROGRAMS = TEOBResumS
TEOBResumS_SOURCES = main.cpp
TEOBResumS_LDADD = libTEOBResumS.la
