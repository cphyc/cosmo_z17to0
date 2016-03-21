CC=ifort
#CFLAGS=-DBIG_RUN -fast -free -m64 -cpp
CFLAGS=-i4 -r8 -O0  -fno-alias -fno-fnalias  -g -debug -traceback -check all  -implicitnone -warn all  -fpe0 -fp-stack-check -ftrapuv -heap-arrays -gen-interface -warn interface -fp-stack-check -ftrapuv
LFLAGS=
SRCFILES=$(shell find . -type f -name "*.f90" -not -path "*tools*")

BINFOLDER=bins
BINFILES=$(patsubst %.f90,$(BINFOLDER)/%,$(SRCFILES))

LINKFILES=$(shell find ./tools -type f -name "*.f90")
MODFILES=$(patsubst %.f90,%.mod,$(LINKFILES))

all: $(MODFILES) $(BINFILES)

%.mod: %.f90 Makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BINFOLDER)/%: %.f90 Makefile $(MODFILES)
	@mkdir -p $(BINFOLDER)
	$(CC) $(CFLAGS) $(MODFILES) -o $@ $<

todolist:
	-@for file in $(SRCFILES:Makefile=); do fgrep -H -e TODO -e FIXME $$file; done; true


clean:
	$(RM) $(BINFILES) $(MODFILES)
