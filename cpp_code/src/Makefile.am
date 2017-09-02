AM_LDFLAGS=`gsl-config --libs` -lgsl -lgslcblas -lm
ACLOCAL_AMFLAGS='-I m4'
bin_PROGRAMS = TEOBResum
TEOBResum_SOURCES = main.cpp

lib_LTLIBRARIES = libTEOBResum.la
libTEOBResum_la_SOURCES =           \
    TEOBResumFlux.cpp               \
    TEOBResumMetric.cpp             \
    TEOBResumUtils.cpp              \
    TEOBResum.cpp                   \
    TEOBResumHlm.cpp                \
    TEOBResumNQC.cpp                \
    TEOBResumWaveform.cpp           \
    TEOBResumInitialConditions.cpp  \
    TEOBResumRingdown.cpp
