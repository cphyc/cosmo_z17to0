F2PY=f2py
CCNAME=intelem
FCNAME=intelem
MODULENAME=_tools
SIGFILE=tools.pyf
SKIPFUN=parse_params read_info_headers read_brick_header read_brick_data read_particle read_mergertree_headers_1 read_mergertree_headers_2 read_mergertree_parent_of read_1 read_1_dynamics read_2_dummy read_2 read_region
SOURCES=compute.f90 misc.f90 io.f90
LIBS=$(shell pwd)/../third_party/FLAP/shared/libflap.so
INCLUDEDIR=../dist-no-p/mod

CFLAGS=-Ofast -fp-model fast=2 -static -free -cpp

all: _tools.so convolution.so readtree.so

_tools.pyf: $(SOURCES)
	$(F2PY) -h $@ --overwrite-signature -m $(basename $@) $^ skip: $(SKIPFUN) :  \
		> /dev/null

_tools.so: _tools.pyf $(SOURCES)
	$(F2PY) -c $^ $(LIBS) -I$(INCLUDEDIR) -m $(basename $@) skip: $(SKIPFUN) : \
		--compiler=$(CCNAME) --fcompiler=$(FCNAME) --f90flags="$(CFLAGS)" \
		> /dev/null

convolution.so: readtree.f90
	$(F2PY) -m $(basename $@) -c $^ \
		--fcompiler=$(FCNAME) --compiler=$(CCNAME) --f90flags="$(CFLAGS)" \
		> /dev/null

readtree.so: readtree.f90
	$(F2PY) -m $(basename $@) -c $^ \
		--fcompiler=$(FCNAME) --compiler=$(CCNAME) --f90flags="$(CFLAGS)" \
		> /dev/null


clean:
	rm -f *.so *.pyc *.pyf
