FORT = gfortran
FLAGS = -Wunused-parameter -mcmodel=medium -O3

version := $(lastword $(shell head -1 illum/__init__.py))

SRCS := $(wildcard illum/kernel/*.f90)
BINS := $(SRCS:illum/kernel/%.f90=bin/%)
LIBS := $(wildcard illum/kernel/libs/*.f90)
F2PY := $(wildcard illum/compute/*.f90)
PYSO := illum/compute/compute.so

.PHONY: all kernel f2py clean

all: kernel f2py


kernel: ${BINS}

bin/illumina: illum/kernel/illumina.f90
	@sed "s/__version__/${version}/" $^ > illumina.f90
	${FORT} ${FLAGS} illumina.f90 ${LIBS} -o $@
	@rm illumina.f90

bin/%: illum/kernel/%.f90
	${FORT} ${FLAGS} $^ -o $@


f2py: ${PYSO}

%.so: ${F2PY}
	f2py3 -c ${F2PY} -m compute
	@mv compute.*.so $@


clean:
	rm -f ${BINS} ${PYSO}
