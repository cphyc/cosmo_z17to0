CC=ifort
CFLAGS=-O3
LFLAGS=
SRCFILES=$(shell find . -type f -name "\*.f90")
BINFILES=$(patsubst %.f90,%,$(SRCFILES))

%: %.f90 Makefile
	$(CC) $(CFLAGS) -o $@ $<

todolist:
	-@for file in $(ALLFILES:Makefile=); do fgrep -H -e TODO -e FIXME $$file; done; true


clean:
	$(RM) $(BINFILES)
