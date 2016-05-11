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

# F2PY=f2py
# CC=icc
# FC=ifort
# CCNAME=intelem
# FCNAME=intelem
# MODULENAME=_tools
# SRCFILES=$(shell find . -type f -name "*.f90" -not -path "./tools/*" -not -path "./f90/*" -not -path "*test*")
# TESTSRCFILES=$(shell find . -type f -path "*test*")
# BINFILES=$(patsubst %.f90,%.so,$(SRCFILES))
# MODFILES=$(patsubst %.f90,%.mod,$(SRCFILES))

# all: $(SRCFILES)
# 	$(F2PY) -c $^ -m $(MODULENAME) --fcompiler=$(FCNAME) --compiler=$(CCNAME)

# modules: $(MODFILES)

# %.mod: %.f90 Makefile
# 	$(FC) -c $< -o $@

# # %.so: %.f90 Makefile
# # 	$(F2PY) -c $< -m $(MODULENAME) --fcompiler=$(FCNAME) --compiler=$(CCNAME)

# test:
# 	$(FC) $(TESTSRCFILES) -o test
# 	./test

# clean:
# 	$(RM) $(BINFILES) $(MODFILES)
