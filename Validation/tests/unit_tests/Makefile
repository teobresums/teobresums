# Makefile for spin dynamics tests

# make sure to use bash on Debian/Ubuntu, not dash
SHELL=/bin/bash

#
# names and dirs 
#

NAME=unit_tests
EXE=$(NAME).x
BASE=$(shell /bin/pwd)
SRC=../../../C/src/TEOBResumSUtils.c\
 ../../../C/src/TEOBResumSPars.c\
 ../../../C/src/TEOBResumSDynamics.c\
 ../../../C/src/TEOBResumSFits.c\
 ../../../C/src/TEOBResumSFlux.c\
 ../../../C/src/TEOBResumSInitialConditions.c\
 ../../../C/src/TEOBResumSMetric.c\
 ../../../C/src/TEOBResumSOMP.c\
 ../../../C/src/TEOBResumSPostAdiab.c\
 ../../../C/src/TEOBResumSWaveform.c\
 ../../../C/src/SpecialFuns.c\
 unit_tests.c

#
# compilation options
#

CC=gcc
CFLAGS=-std=gnu99 `gsl-config --cflags` # Debug
#CFLAGS=-std=c99 `gsl-config --cflags` # Debug

CFLAGS +=-g  # Debug
#CFLAGS +=-pg # Debug & Profile
#CFLAGS +=-O3  # Optimized

CFLAGS += -DDEBUG=0 # 1/0 to turn on/off debug mode
CFLAGS += -DVERBOSE=0 # 1/0 to turn on/off verbose mode

#CFLAGS += -fopenmp # OpenMP
#CFLAGS += -DUSETIMERS=0 # OpenMP timers, note they need OpenMP

CFLAGS += -DUSEGRAVITOMAGNETICTERMS=1 # 1/0 to turn on gravitomagnetic terms
CFLAGS += -DUSEBTIDALPOTENTIAL=1 # 1/0 to turn on B tidal terms
#CFLAGS += -DUSERK45=0 # use rk45? rk8pd is default
#CFLAGS += -DPOSTPOSTCIRCULAR=0 # 1/0 to turn on post-post-circular initial conditions
#CFLAGS += -DEXCLUDESPINSPINTIDES=1 # 1/0 to turn off/on (exclude) spin^2-tidal couplings
#CFLAGS += -DEOBRUN=1 # use old/new main #TO BE REMOVED (use only new)

LFLAGS=`gsl-config --libs`
LDLFLAGS=-lgsl -lgslcblas -lm -lconfig

#
# rules
#

all: $(EXE) 
	@echo "All done"

$(EXE): $(SRC)
	@echo "Building $@ ..."
	@$(CC) $(CFLAGS) $(LFLAGS) $(SRC) -o $@ $(LDLFLAGS)

new: clean $(EXE)
	@echo "All done"

.PHONY: all new 

# DO NOT DELETE

