# makefile to compile and run benchmark

# Add a parfile in the list below to run in in the benchmark
# All these parfile are and should be kept in ../C/par/
FILES=bbh_q1_s0s0_10Hz.par bbh_q1_s0s0_10Hz_postadiab.par bbh_q1_s0s0_M10_10Hz_postadiab.par bbh_q3_s08s-02_15Hz_postadiab.par bbh_q3_s08s-02_M12_15Hz.par bbh_q3_s08s-02_M12_15Hz_postadiab.par bns_q1_s01s01_20Hz.par bns_q1_s01s01_20Hz_postadiab.par bns_q1_s0s0_10Hz.par bns_q1_s0s0_10Hz_postadiab.par bns_q1_s0s0_M2.7_10Hz.par 

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

