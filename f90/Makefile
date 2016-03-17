F2PY=f2py
CC=icc
FC=ifort
CCNAME=intelem
FCNAME=intelem
MODULENAME=tools
SRCFILES=$(shell find . -type f -name "*.f90" -not -path "./tools/*" -not -path "./f90/*")
BINFILES=$(patsubst %.f90,%.so,$(SRCFILES)) $(patsubst %.f90,%.o,$(SRCFILES))

all: $(BINFILES)

%.so: %.f90 Makefile
	$(F2PY) -c $< -m $(MODULENAME) --fcompiler=$(FCNAME) --compiler=$(CCNAME)

clean:
	$(RM) $(BINFILES) $(MODFILES)
