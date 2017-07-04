CC=g++ -O4
CFLAGS=-c -std=c++11
INC=-I/opt/local/include
LIBS=-L/opt/local/lib
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
input_struc.o 	     \
interp_grid.o 	     \
interpolate_wf.o     \
interpolator.o	     \
main.o

all:main

main:$(OBJS)
	$(CC) $(LIBS)  $(OBJS) -o EOB_ihes.out $(LDLIBS)
main.o: main.cpp
	$(CC) $(CFLAGS) $(INC) main.cpp
A_NumDenom.o: A_NumDenom.cpp
	$(CC) $(CFLAGS) $(INC) A_NumDenom.cpp
AdiabLR.o: AdiabLR.cpp
	$(CC) $(CFLAGS) $(INC) AdiabLR.cpp
deltalm.o: deltalm.cpp
	$(CC) $(CFLAGS) $(INC) deltalm.cpp
dtnqc_fit.o: dtnqc_fit.cpp
	$(CC) $(CFLAGS) $(INC) dtnqc_fit.cpp
Eulerlog.o: Eulerlog.cpp
	$(CC) $(CFLAGS) $(INC) Eulerlog.cpp
f_lm.o: f_lm.cpp
	$(CC) $(CFLAGS) $(INC) f_lm.cpp
FDdrvt_omega.o: FDdrvt_omega.cpp
	$(CC) $(CFLAGS) $(INC) FDdrvt_omega.cpp
FDdrvt.o: FDdrvt.cpp
	$(CC) $(CFLAGS) $(INC) FDdrvt.cpp
file_names.o: file_names.cpp
	$(CC) $(CFLAGS) $(INC) file_names.cpp
find_a1a2a3.o: find_a1a2a3.cpp
	$(CC) $(CFLAGS) $(INC) find_a1a2a3.cpp
FlmNewt.o: FlmNewt.cpp
	$(CC) $(CFLAGS) $(INC) FlmNewt.cpp
flux.o: flux.cpp
	$(CC) $(CFLAGS) $(INC) flux.cpp
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
hlmNQC.o: hlmNQC.cpp
	$(CC) $(CFLAGS) $(INC) hlmNQC.cpp
HorizonFlux.o: HorizonFlux.cpp
	$(CC) $(CFLAGS) $(INC) HorizonFlux.cpp
initial.o: initial.cpp
	$(CC) $(CFLAGS) $(INC) initial.cpp
input_struc.o: input_struc.cpp
	$(CC) $(CFLAGS) $(INC) input_struc.cpp
interp_grid.o: interp_grid.cpp
	$(CC) $(CFLAGS) $(INC) interp_grid.cpp
interpolate_wf.o: interpolate_wf.cpp
	$(CC) $(CFLAGS) $(INC) interpolate_wf.cpp
interpolator.o: interpolator.cpp
	$(CC) $(CFLAGS) $(INC) interpolator.cpp

clean:
	rm -rf *.o EOB_ihes.out
