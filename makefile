FORT = gfortran
FLAGS = -Wunused-parameter -mcmodel=medium -O3

version := $(lastword $(shell head -1 illum/__init__.py))

SRCS := $(wildcard illum/kernel/*.f)
BINS := $(SRCS:illum/kernel/%.f=bin/%)
LIBS := $(wildcard illum/kernel/libs/*.f)
F2PY := $(wildcard illum/compute/*.f90)
PYSO := illum/compute/compute.so

.PHONY: all f2py clean

all: kernel f2py


kernel: ${BINS}

bin/illumina: illum/kernel/illumina.f
	@sed "s/__version__/${version}/" $^ > illumina.f
	${FORT} ${FLAGS} illumina.f ${LIBS} -o $@
	@rm illumina.f

bin/%: illum/kernel/%.f
	${FORT} ${FLAGS} $^ -o $@


f2py: ${PYSO}

%.so: ${F2PY}
	f2py3 -c -DF2PY_REPORT_ON_ARRAY_COPY=1 ${F2PY} -m compute
	@mv compute.*.so $@


clean:
	rm -f ${BINS} ${PYSO}
