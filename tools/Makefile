F2PY=f2py
CCNAME=intelem
FCNAME=intelem
MODULENAME=_tools
SIGFILE=tools.pyf
SKIPFUN=parse_params read_info_headers read_brick_header read_brick_data read_particle read_mergertree_headers_1 read_mergertree_headers_2 read_mergertree_parent_of read_1 read_1_dynamics read_2_dummy read_2
SOURCES=misc.f90 io.f90
LIBS=$(shell pwd)/../third_party/FLAP/shared/libflap.so
INCLUDEDIR=../dist-no-p/mod

all: compile

tools.pyf: $(SOURCES)
	$(F2PY) -h $@ --overwrite-signature -m $(MODULENAME) misc.f90 io.f90 skip: $(SKIPFUN) :

compile: tools.pyf $(SOURCES)
	$(F2PY) -c $^ $(LIBS) -I$(INCLUDEDIR) --compiler=$(CCNAME) --fcompiler=$(FCNAME)
