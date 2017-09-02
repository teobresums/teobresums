CC=gcc -g
CFLAGS=-c -Wall -std=c99
INC=`gsl-config --cflags`
LIBS=`gsl-config --libs`
LDLIBS=-lgsl -lgslcblas -lm

OBJS =               \
A_NumDenom.o	     \
AdiabLR.o    	     \
deltalm.o    	     \
dtnqc_fit.o 	     \
Eulerlog.o  	     \
f_lm.o	     	     \
FDdrvt_omega.o 	     \
FDdrvt.o 	     \
file_names.o	     \
find_a1a2a3.o 	     \
FlmNewt.o 	     \
flux.o 		     \
get_Omg_orb.o 	     \
HealyBBHFitRemnant.o \
hhatlmtail.o 	     \
hlm_Tidal.o 	     \
hlm.o 		     \
hlmNewt.o 	     \
hlmNQC.o 	     \
HorizonFlux.o	     \
initial.o 	     \
interp_grid.o 	     \
interpolator.o	     \
main.o         \
Metric.o       \
NQC.o         \
QNMHybridFitCab.o  \
read_config.o	     \
RHS.o           \
s_A5PNlog.o	     \
s_bisec.o        \
s_D1.o		     \
s_flm.o		     \
s_Flux.o	     \
s_get_rc.o	     \
s_GS.o		     \
s_HorizonFlux.o	     \
s_initial.o	     \
s_Metric.o	     \
s_RHS.o		     \
s_waveform.o	    \
SpeedyTail.o      \
Tlm.o     \
constants.o \
multipole_index.o  \
Array.o

#interpolate_wf.o     \
#Q_omega.o     \
#ringdown.o	     \
#ringdown_match.o     \



all:main

main:$(OBJS)
	$(CC) $(LIBS)  $(OBJS) -o EOB_ihes.out $(LDLIBS)
main.o: main.c
	$(CC) $(CFLAGS) $(INC) main.c
A_NumDenom.o: A_NumDenom.c
	$(CC) $(CFLAGS) $(INC) A_NumDenom.c
AdiabLR.o: AdiabLR.c
	$(CC) $(CFLAGS) $(INC) AdiabLR.c
deltalm.o: deltalm.c
	$(CC) $(CFLAGS) $(INC) deltalm.c
dtnqc_fit.o: dtnqc_fit.c
	$(CC) $(CFLAGS) $(INC) dtnqc_fit.c
Eulerlog.o: Eulerlog.c
	$(CC) $(CFLAGS) $(INC) Eulerlog.c
f_lm.o: f_lm.c
	$(CC) $(CFLAGS) $(INC) f_lm.c
FDdrvt_omega.o: FDdrvt_omega.c
	$(CC) $(CFLAGS) $(INC) FDdrvt_omega.c
FDdrvt.o: FDdrvt.c
	$(CC) $(CFLAGS) $(INC) FDdrvt.c
file_names.o: file_names.c
	$(CC) $(CFLAGS) $(INC) file_names.c
find_a1a2a3.o: find_a1a2a3.c
	$(CC) $(CFLAGS) $(INC) find_a1a2a3.c
FlmNewt.o: FlmNewt.c
	$(CC) $(CFLAGS) $(INC) FlmNewt.c
flux.o: flux.c
	$(CC) $(CFLAGS) $(INC) flux.c
get_Omg_orb.o: get_Omg_orb.c
	$(CC) $(CFLAGS) $(INC) get_Omg_orb.c
HealyBBHFitRemnant.o: HealyBBHFitRemnant.c
	$(CC) $(CFLAGS) $(INC) HealyBBHFitRemnant.c
hhatlmtail.o: hhatlmtail.c
	$(CC) $(CFLAGS) $(INC) hhatlmtail.c
hlm_Tidal.o: hlm_Tidal.c
	$(CC) $(CFLAGS) $(INC) hlm_Tidal.c
hlm.o: hlm.c
	$(CC) $(CFLAGS) $(INC) hlm.c
hlmNewt.o: hlmNewt.c
	$(CC) $(CFLAGS) $(INC) hlmNewt.c
hlmNQC.o: hlmNQC.c
	$(CC) $(CFLAGS) $(INC) hlmNQC.c
HorizonFlux.o: HorizonFlux.c
	$(CC) $(CFLAGS) $(INC) HorizonFlux.c
initial.o: initial.c
	$(CC) $(CFLAGS) $(INC) initial.c
interp_grid.o: interp_grid.c
	$(CC) $(CFLAGS) $(INC) interp_grid.c
interpolate_wf.o: interpolate_wf.c
	$(CC) $(CFLAGS) $(INC) interpolate_wf.c
interpolator.o: interpolator.c
	$(CC) $(CFLAGS) $(INC) interpolator.c
Metric.o: Metric.c
	$(CC) $(CFLAGS) $(INC) Metric.c
NQC.o: NQC.c
	$(CC) $(CFLAGS) $(INC) NQC.c
Q_omega.o: Q_omega.c
	$(CC) $(CFLAGS) $(INC) Q_omega.c
QNMHybridFitCab.o: QNMHybridFitCab.c
	$(CC) $(CFLAGS) $(INC) QNMHybridFitCab.c
read_config.o: read_config.c
	$(CC) $(CFLAGS) $(INC) read_config.c
RHS.o: RHS.c
	$(CC) $(CFLAGS) $(INC) RHS.c
ringdown.o: ringdown.c
	$(CC) $(CFLAGS) $(INC) ringdown.c
ringdown_match.o: ringdown_match.c
	$(CC) $(CFLAGS) $(INC) ringdown_match.c
s_A5PNlog.o: s_A5PNlog.c
	$(CC) $(CFLAGS) $(INC) s_A5PNlog.c
s_bisec.o: s_bisec.c
	$(CC) $(CFLAGS) $(INC) s_bisec.c
s_D1.o: s_D1.c
	$(CC) $(CFLAGS) $(INC) s_D1.c
s_flm.o: s_flm.c
	$(CC) $(CFLAGS) $(INC) s_flm.c
s_Flux.o: s_Flux.c
	$(CC) $(CFLAGS) $(INC) s_Flux.c
s_get_rc.o: s_get_rc.c
	$(CC) $(CFLAGS) $(INC) s_get_rc.c
s_GS.o: s_GS.c
	$(CC) $(CFLAGS) $(INC) s_GS.c
s_HorizonFlux.o: s_HorizonFlux.c
	$(CC) $(CFLAGS) $(INC) s_HorizonFlux.c
s_initial.o: s_initial.c
	$(CC) $(CFLAGS) $(INC) s_initial.c
s_Metric.o: s_Metric.c
	$(CC) $(CFLAGS) $(INC) s_Metric.c
s_RHS.o: s_RHS.c
	$(CC) $(CFLAGS) $(INC) s_RHS.c
s_waveform.o: s_waveform.c
	$(CC) $(CFLAGS) $(INC) s_waveform.c
SpeedyTail.o: SpeedyTail.c
	$(CC) $(CFLAGS) $(INC) SpeedyTail.c
Tlm.o: Tlm.c
	$(CC) $(CFLAGS) $(INC) Tlm.c

constants.o: constants.c
	$(CC) $(CFLAGS) $(INC) constants.c
multipole_index.o: multipole_index.c
	$(CC) $(CFLAGS) $(INC) multipole_index.c
Array.o: Array.c
	$(CC) $(CFLAGS) $(INC) Array.c



clean:
	rm -rf *.o EOB_ihes.out

