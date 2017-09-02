CC=g++
CFLAGS=-c -std=c++11 -g -Wall -fPIC
INC=`gsl-config --cflags`
LIBS=`gsl-config --libs`
LDLIBS=-lgsl -lgslcblas -lm

OBJS =               \
main.o		\
TEOBResum.o  \
TEOBResumFlux.o \
TEOBResumHlm.o \
TEOBResumInitialConditions.o \
TEOBResumMetric.o \
TEOBResumNQC.o	\
TEOBResumRingdown.o \
TEOBResumUtils.o \
TEOBResumWaveform.o \

all:main

main:$(OBJS)
	$(CC) $(LIBS)  $(OBJS) -o EOB_ihes.out $(LDLIBS)
main.o: main.cpp
	$(CC) $(CFLAGS) $(INC) main.cpp
TEOBResum.o: TEOBResum.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResum.cpp
TEOBResumFlux.o: TEOBResumFlux.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumFlux.cpp
TEOBResumHlm.o: TEOBResumHlm.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumHlm.cpp
TEOBResumInitialConditions.o: TEOBResumInitialConditions.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumInitialConditions.cpp
TEOBResumMetric.o: TEOBResumMetric.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumMetric.cpp
TEOBResumNQC.o: TEOBResumNQC.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumNQC.cpp
TEOBResumRingdown.o: TEOBResumRingdown.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumRingdown.cpp
TEOBResumUtils.o: TEOBResumUtils.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumUtils.cpp
TEOBResumWaveform.o: TEOBResumWaveform.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumWaveform.cpp

clean:
	rm -rf *.o *.h.gch EOB_ihes.out

