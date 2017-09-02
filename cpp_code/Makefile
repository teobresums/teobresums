CC=g++
CFLAGS=-c -std=c++11 -g -Wall -fPIC
INC=`gsl-config --cflags`
LIBS=`gsl-config --libs`
LDLIBS=-lgsl -lgslcblas -lm

OBJS =               \
A_NumDenom.o	     \
get_Omg_orb.o 	     \
LALSimIMRTEOBIHES.o  \
main.o		\
read_config.o	     \
s_bisec.o        \
s_D1.o		     \
s_get_rc.o	     \
s_GS.o		     \
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
LALSimIMRTEOBIHES.o: LALSimIMRTEOBIHES.cpp
	$(CC) $(CFLAGS) $(INC) LALSimIMRTEOBIHES.cpp
A_NumDenom.o: A_NumDenom.cpp
	$(CC) $(CFLAGS) $(INC) A_NumDenom.cpp
get_Omg_orb.o: get_Omg_orb.cpp
	$(CC) $(CFLAGS) $(INC) get_Omg_orb.cpp
TEOBResumHlm.o: TEOBResumHlm.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumHlm.cpp
read_config.o: read_config.cpp
	$(CC) $(CFLAGS) $(INC) read_config.cpp
s_bisec.o: s_bisec.cpp
	$(CC) $(CFLAGS) $(INC) s_bisec.cpp
s_D1.o: s_D1.cpp
	$(CC) $(CFLAGS) $(INC) s_D1.cpp
s_get_rc.o: s_get_rc.cpp
	$(CC) $(CFLAGS) $(INC) s_get_rc.cpp
s_GS.o: s_GS.cpp
	$(CC) $(CFLAGS) $(INC) s_GS.cpp
TEOBResumFlux.o: TEOBResumFlux.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumFlux.cpp
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

