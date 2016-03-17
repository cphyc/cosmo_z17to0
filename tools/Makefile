F2PY=f2py
CC=icc
FC=ifort
CCNAME=intelem
FCNAME=intelem
MODULENAME=_tools
SRCFILES=$(shell find . -type f -name "*.f90" -not -path "./tools/*" -not -path "./f90/*" -not -path "*test*")
TESTSRCFILES=$(shell find . -type f -path "*test*")
BINFILES=$(patsubst %.f90,%.so,$(SRCFILES))

all: $(BINFILES)

%.so: %.f90 Makefile
	$(F2PY) -c $< -m $(MODULENAME) --fcompiler=$(FCNAME) --compiler=$(CCNAME)

test:
	$(FC) $(TESTSRCFILES) -o test
	./test

clean:
	$(RM) $(BINFILES) $(MODFILES)
