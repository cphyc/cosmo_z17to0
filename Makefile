CC=ifort
CFLAGS=
LFLAGS=

%:
	$(CC) -o $* $*.f90
