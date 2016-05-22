.SUFFIXES:
.SUFFIXES:  .out .o .f90 .for .f
# Edit this to be correct on your system
F90 = ifort
MODULEDIR= include

ifeq ($(F90),ifort)
	F90OPTS = -O3 -mcmodel=medium -module $(MODULEDIR) -I$(MODULEDIR) -qopenmp -fPIC
	LDOPTS =

	MKLARCH = em64t
	MKLPATH=$(MKLROOT)/lib/$(MKLARCH)

	LAPACK  = -I$(MKLROOT)/include -mkl

	LAPACKLIB = $(LAPACK) -lpthread
endif
#-----------------------------------------------------

OBJ = extrema_types.o extrema_storage.o extrema_mod.o

objects: $(OBJ)
#----------------------------------------------------

%.o: %.f90
	$(F90) $(F90OPTS) -c $^

#----------------------------------------------------------

clean:
	/bin/rm -f *.o

cleanall:
	/bin/rm -f *.o $(MODULEDIR)/[!n]*.mod
