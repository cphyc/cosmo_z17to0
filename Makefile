CC=ifort
CFLAGS=-traceback -g
LFLAGS=
SRCFILES=$(shell find . -type f -name "*.f90" -not -path "./tools/*")

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
