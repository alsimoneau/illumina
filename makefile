FORT = gfortran

version := $(lastword $(shell head -1 illum/__init__.py))

F2PY := $(wildcard illum/compute/*.f90)
PYSO := illum/compute/compute.so
LIBS := $(wildcard illum/compute/libs/*.f90)
OBJS := $(LIBS:illum/compute/libs/%.f90=lib/%.o)
DEP_FILE := lib/depends.mk

.PHONY: all depend libs f2py clean

all: depend libs f2py

depend: $(DEP_FILE)

$(DEP_FILE): $(LIBS)
	@echo "Generating dependencies tree..."
	@mkdir -p $(@D)
	@python bin/fortdepends.py $(LIBS)

include $(DEP_FILE)

libs: $(OBJS)

lib/%.o: illum/compute/libs/%.f90
	$(FORT) -fPIC -Jlib -O3 -c $< -o $@

f2py: ${PYSO}

%.so: ${F2PY}
	f2py3 --fcompiler=$(FORT) --f90flags='-fopenmp' -lgomp -Ilib -c ${F2PY} ${OBJS} -m compute
	@mv compute.*.so $@

clean:
	rm -rf lib
	rm -f ${PYSO}
