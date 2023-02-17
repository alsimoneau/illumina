FORT = gfortran

version := $(lastword $(shell head -1 illum/__init__.py))

F2PY := $(wildcard illum/compute/*.f90)
PYSO := illum/compute/compute.so
LIBS := $(wildcard illum/compute/libs/*.f90)
OBJS := $(LIBS:illum/compute/libs/%.f90=lib/%.o)

.PHONY: all libs f2py clean

all: libs f2py

f2py: ${PYSO}

libs: lib/math.o $(OBJS)

lib/%.o: illum/compute/libs/%.f90
	@mkdir -p $(@D)
	$(FORT) -fPIC -Jlib -O3 -c $< -o $@

%.so: ${F2PY}
	f2py3 --fcompiler=$(FORT) --f90flags='-fopenmp' -lgomp -Ilib -c ${F2PY} ${OBJS} -m compute
	@mv compute.*.so $@

print:  
	@echo $(MODS)

clean:
	rm -rf lib
	rm -f ${PYSO}
