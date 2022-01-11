# makefile to compile and run benchmark
# uses parfiles in Python/Utils/parfiles.json

# do not touch below here

TOP=$(shell pwd)
SRC=$(abspath ../../C/)
UT=$(abspath ../)
MAKEFILE=$(SRC)/Makefile
NOW:=$(shell date +'%y%m%d')
GIT:=$(shell git rev-list HEAD | wc -l)

WORKDIR=$(TOP)/bench_$(GIT)/$(NOW)

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
	export TEOBRESUMS=$(WORKDIR)&&\
	cd $(WORKDIR) &&\
	cp $(UT)/eobutils.py . && cp $(UT)/parfile.py . && cp $(UT)/parfiles.json . &&\
	cp $(TOP)/benchmark.py . &&\
	python3 benchmark.py &&\
	rm *.py

clean: 
	cd $(WORKDIR) && rm *.x

