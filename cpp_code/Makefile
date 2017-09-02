CC=g++
CFLAGS=-c -std=c++11 -g -Wall -fPIC
INC=`gsl-config --cflags`
LIBS=`gsl-config --libs`
LDLIBS=-lgsl -lgslcblas -lm

OBJS =               \
A_NumDenom.o	     \
AdiabLR.o    	     \
deltalm.o    	     \
f_lm.o	     	     \
FlmNewt.o 	     \
get_Omg_orb.o 	     \
HealyBBHFitRemnant.o \
hhatlmtail.o 	     \
hlm_Tidal.o 	     \
hlm.o 		     \
hlmNewt.o 	     \
HorizonFlux.o	     \
LALSimIMRTEOBIHES.o  \
main.o		\
Metric.o       \
read_config.o	     \
RHS.o           \
s_A5PNlog.o	     \
s_bisec.o        \
s_D1.o		     \
s_flm.o		     \
s_get_rc.o	     \
s_GS.o		     \
s_HorizonFlux.o	     \
s_Metric.o	     \
s_RHS.o		     \
s_waveform.o	    \
TEOBResumFlux.o \
TEOBResumInitialConditions.o \
TEOBResumNQC.o	\
TEOBResumRingdown.o \
TEOBResumUtils.o \
Tlm.o \

all:main

main:$(OBJS)
	$(CC) $(LIBS)  $(OBJS) -o EOB_ihes.out $(LDLIBS)
main.o: main.cpp
	$(CC) $(CFLAGS) $(INC) main.cpp
LALSimIMRTEOBIHES.o: LALSimIMRTEOBIHES.cpp
	$(CC) $(CFLAGS) $(INC) LALSimIMRTEOBIHES.cpp
A_NumDenom.o: A_NumDenom.cpp
	$(CC) $(CFLAGS) $(INC) A_NumDenom.cpp
AdiabLR.o: AdiabLR.cpp
	$(CC) $(CFLAGS) $(INC) AdiabLR.cpp
deltalm.o: deltalm.cpp
	$(CC) $(CFLAGS) $(INC) deltalm.cpp
f_lm.o: f_lm.cpp
	$(CC) $(CFLAGS) $(INC) f_lm.cpp
FlmNewt.o: FlmNewt.cpp
	$(CC) $(CFLAGS) $(INC) FlmNewt.cpp
get_Omg_orb.o: get_Omg_orb.cpp
	$(CC) $(CFLAGS) $(INC) get_Omg_orb.cpp
HealyBBHFitRemnant.o: HealyBBHFitRemnant.cpp
	$(CC) $(CFLAGS) $(INC) HealyBBHFitRemnant.cpp
hhatlmtail.o: hhatlmtail.cpp
	$(CC) $(CFLAGS) $(INC) hhatlmtail.cpp
hlm_Tidal.o: hlm_Tidal.cpp
	$(CC) $(CFLAGS) $(INC) hlm_Tidal.cpp
hlm.o: hlm.cpp
	$(CC) $(CFLAGS) $(INC) hlm.cpp
hlmNewt.o: hlmNewt.cpp
	$(CC) $(CFLAGS) $(INC) hlmNewt.cpp
HorizonFlux.o: HorizonFlux.cpp
	$(CC) $(CFLAGS) $(INC) HorizonFlux.cpp
Metric.o: Metric.cpp
	$(CC) $(CFLAGS) $(INC) Metric.cpp
read_config.o: read_config.cpp
	$(CC) $(CFLAGS) $(INC) read_config.cpp
RHS.o: RHS.cpp
	$(CC) $(CFLAGS) $(INC) RHS.cpp
s_A5PNlog.o: s_A5PNlog.cpp
	$(CC) $(CFLAGS) $(INC) s_A5PNlog.cpp
s_bisec.o: s_bisec.cpp
	$(CC) $(CFLAGS) $(INC) s_bisec.cpp
s_D1.o: s_D1.cpp
	$(CC) $(CFLAGS) $(INC) s_D1.cpp
s_flm.o: s_flm.cpp
	$(CC) $(CFLAGS) $(INC) s_flm.cpp
s_get_rc.o: s_get_rc.cpp
	$(CC) $(CFLAGS) $(INC) s_get_rc.cpp
s_GS.o: s_GS.cpp
	$(CC) $(CFLAGS) $(INC) s_GS.cpp
s_HorizonFlux.o: s_HorizonFlux.cpp
	$(CC) $(CFLAGS) $(INC) s_HorizonFlux.cpp
s_Metric.o: s_Metric.cpp
	$(CC) $(CFLAGS) $(INC) s_Metric.cpp
s_RHS.o: s_RHS.cpp
	$(CC) $(CFLAGS) $(INC) s_RHS.cpp
s_waveform.o: s_waveform.cpp
	$(CC) $(CFLAGS) $(INC) s_waveform.cpp
spinsphericalharm.o: spinsphericalharm.cpp
	$(CC) $(CFLAGS) $(INC) spinsphericalharm.cpp
TEOBResumFlux.o: TEOBResumFlux.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumFlux.cpp
TEOBResumInitialConditions.o: TEOBResumInitialConditions.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumInitialConditions.cpp
TEOBResumNQC.o: TEOBResumNQC.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumNQC.cpp
TEOBResumRingdown.o: TEOBResumRingdown.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumRingdown.cpp
TEOBResumUtils.o: TEOBResumUtils.cpp
	$(CC) $(CFLAGS) $(INC) TEOBResumUtils.cpp
Tlm.o: Tlm.cpp
	$(CC) $(CFLAGS) $(INC) Tlm.cpp

clean:
	rm -rf *.o *.h.gch EOB_ihes.out

