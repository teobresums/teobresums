# makefile to compile and run benchmark
TOP=$(shell pwd)
SRC=$(abspath ../C++/)
RUN=$(abspath ./runtest.sh)
PAR=$(abspath ./parfiles)
MAKEFILE=$(SRC)/Makefile.TEOBResumS
NOW:=$(shell date +'%y%m%d')
GIT:=$(shell git rev-list HEAD | wc -l)

WORKDIR=$(TOP)/test_$(GIT)/$(NOW)

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
	cd $(WORKDIR) && \
	cp -r $(PAR) . && \
	cp $(RUN) . &&\
	./runtest.sh

clean: 
	cd $(WORKDIR) && rm *.x
