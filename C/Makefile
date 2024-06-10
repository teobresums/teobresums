# Makefile for TEOBResumS

# make sure to use bash on Debian/Ubuntu, not dash
SHELL=/bin/bash

#
# names and dirs 
#

BASE=$(shell /bin/pwd)
SRCD=$(BASE)/src
OBJD=$(BASE)/obj
LIBD=$(BASE)/lib

NAME=TEOBResumS
EXE=$(NAME).x
SRC=$(wildcard $(SRCD)/*.c)
OBJ=$(patsubst $(SRCD)/%.c,$(OBJD)/%.o,$(SRC))
LIB=$(LIBD)/lib$(NAME).a

#
# compilation options
#

CC=gcc
CFLAGS=-std=c99 `gsl-config --cflags` # Debug

CFLAGS +=-g  # Debug
#CFLAGS +=-pg # Debug & Profile
#CFLAGS +=-O3  # Optimized

CFLAGS += -DDEBUG=0 # 1/0 to turn on/off debug mode
CFLAGS += -DVERBOSE=0 # 1/0 to turn on/off verbose mode

#CFLAGS += -fcommon # gcc 10
#CFLAGS += -fopenmp # OpenMP
#CFLAGS += -DUSETIMERS=0 # OpenMP timers, note they need OpenMP

CFLAGS += -DUSEBTIDALPOTENTIAL=1 # 1/0 to turn on B tidal terms
#CFLAGS += -DUSERK45=0 # use rk45? rk8pd is default
#CFLAGS += -DPOSTPOSTCIRCULAR=0 # 1/0 to turn on post-post-circular initial conditions
#CFLAGS += -DEXCLUDESPINSPINTIDES=1 # 1/0 to turn off/on (exclude) spin^2-tidal couplings
#CFLAGS += -DUSEGRAVITOMAGNETICTERMS=1 # 1/0 to turn on gravitomagnetic terms

LFLAGS=`gsl-config --libs`
LDLFLAGS=-lgsl -lgslcblas -lm 

#
# rules
#

all: $(EXE) # pywrap $(LIB)
	@echo "All done"
	
$(EXE): $(OBJ)
	@echo "Building $@ ..."
	@$(CC) $(CFLAGS) $(LFLAGS) $(OBJ) -o $@ $(LDLFLAGS)

$(OBJD)/%.o: $(SRCD)/%.c
	@mkdir -p $(dir $@)
	@echo "Compiling $< ... -> $@"
	@$(CC) $(CFLAGS) -c $< -o $@

clean: 
	@echo "Cleaning ..."
	@rm -vf $(OBJ) 
	@rm -vf $(LIB) 
	@echo "... done"

new: clean $(EXE)
	@echo "All done"

clearup: clean 
	@echo "Removing dirs and exe ..."
	@rm -rvf $(OBJ) 
	@rm -rvf $(LIB) 
	@rm -rvf $(EXE) 
	@echo "All done"

.PHONY: all clean new #help

# DO NOT DELETE
