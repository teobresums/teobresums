# makefile to compile and run benchmark

# Add a parfile in the list below to run in in the benchmark
# All these parfile are and should be kept in ../C/par/
FILES=GW150914.par GW170817.par bbh_q1_s0s0_10Hz.par bbh_q3_s-099s-099_10Hz.par bbh_q5_s09s-04_HM.par bhns_q3_s075s0_f0004.par bns_q1_s-01s-01_10Hz.par

# do not touch below here

TOP=$(shell pwd)
SRC=$(abspath ../C/)
RUN=$(abspath ./Run.sh)
PAR=$(abspath ../C/par/)
MAKEFILE=$(SRC)/Makefile.TEOBResumS
NOW:=$(shell date +'%y%m%d')
GIT:=$(shell git rev-list HEAD | wc -l)

WORKDIR=$(TOP)/bench_$(GIT)/$(NOW)
PARFILES=$(patsubst %,$(PAR)/%,$(FILES))

all: build run

build: $(WORKDIR)

$(WORKDIR):
	mkdir -pv $(WORKDIR)
	cd $(WORKDIR) && \
	cp -r $(SRC)/src/ . && cp $(MAKEFILE) Makefile
	$(MAKE) -C $(WORKDIR) all
	cd $(WORKDIR) && \
	tar zcvf src.tgz Makefile src/ && \
	rm -rvf src/ obj/ Makefile && \
	git log --format="%H" -n 1 > git_status.txt 

run: $(WORKDIR)
	export TEOBRESUMS=$(SRC)&&\
	echo $$TEOBRESUM &&\
	cd $(WORKDIR) &&\
	cp $(PARFILES) . &&\
	cp $(RUN) . &&\
	./Run.sh

clean: 
	cd $(WORKDIR) && rm *.x

