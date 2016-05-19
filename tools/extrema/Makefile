.SUFFIXES:
.SUFFIXES:  .out .o .f90 .for .f
# Edit this to be correct on your system
F90 = ifort
MODULEDIR= include

ifeq ($(F90),ifort)
	F90OPTS = -O3 -mcmodel=medium -module $(MODULEDIR) -I$(MODULEDIR) -qopenmp
	LDOPTS =

	MKLARCH = em64t
	MKLPATH=$(MKLROOT)/lib/$(MKLARCH)

	LAPACK  = -I$(MKLROOT)/include -mkl

	LAPACKLIB = $(LAPACK) -lpthread
endif
#-----------------------------------------------------

OBJ = extrema_types.o extrema_mod.o

objects: $(OBJ)
#----------------------------------------------------

extrema_types.o: extrema_types.f90
	$(F90) $(F90OPTS) -c extrema_types.f90
extrema_mod.o:	extrema_mod.f90
	$(F90) $(F90OPTS) -c extrema_mod.f90

#----------------------------------------------------------

clean:
	/bin/rm -f *.o

cleanall:
	/bin/rm -f *.o $(MODULEDIR)/[!n]*.mod
