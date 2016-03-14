CC=ifort
CFLAGS=-O3
LFLAGS=
SRCFILES=$(shell find . -type f -name "*.f90")
BINFOLDER=bins
BINFILES=$(patsubst %.f90,$(BINFOLDER)/%,$(SRCFILES))

all: $(BINFILES)
	echo $(WILDCARD *.f90)

$(BINFOLDER)/%: %.f90 Makefile
	$(CC) $(CFLAGS) -o $@ $<

todolist:
	-@for file in $(ALLFILES:Makefile=); do fgrep -H -e TODO -e FIXME $$file; done; true


clean:
	$(RM) $(BINFILES)
