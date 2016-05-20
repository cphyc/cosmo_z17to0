
F2PY=f2py
CCNAME=intelem
FCNAME=intelem
FCC=ifort
MODULENAME=_tools
SIGFILE=tools.pyf
SKIPFUN=parse_params read_info_headers read_brick_header read_brick_data read_particle read_mergertree_headers_1 read_mergertree_headers_2 read_mergertree_parent_of read_1 read_1_dynamics read_2_dummy read_2 read_region
SOURCES=
LIBS=$(shell pwd)/../third_party/FLAP/shared/libflap.so
INCLUDEDIR=../dist-no-p/mod

MKLARCH = em64t
MKLPATH=$(MKLROOT)/lib/$(MKLARCH)
# Other choices:  32 64
#

FFTWFLAGS=-lfftw3 -L/softs/fftw/3.3.4-pureintel/lib -I/softs/fftw/3.3.4-pureintel/include
MKLFLAGS= -L/softs/fftw/3.3.4-pureintel/lib -I/softs/intel/compilers_and_libraries_2016.3.210/linux/mkl/include
LAPACK  = -I$(MKLROOT)/include -mkl

LAPACKLIB = $(LAPACK) -lpthread

CFLAGS=-Ofast -fast -fp-model fast=2 -static -free -cpp -xHost $(LAPACKLIB) -traceback -mkl -lpthread -I./modules -fPIC

all: _tools.so convolution.so readtree.so extrema.so

_tools.pyf: types.f90 compute.f90 misc.f90 io.f90
	$(F2PY) -h $@ --overwrite-signature -m $(basename $@) $^ skip: $(SKIPFUN) :  \
		> /dev/null

_tools.so: _tools.pyf types.o compute.o misc.o io.o
	$(F2PY) -c $^ $(LIBS) -I$(INCLUDEDIR) -m $(basename $@) skip: $(SKIPFUN) : \
		--compiler=$(CCNAME) --fcompiler=$(FCNAME) --f90flags="$(CFLAGS)" \
		> /dev/null

convolution.so: types.o convolution.f90
	$(F2PY) -m $(basename $@) -c $^ \
		--fcompiler=$(FCNAME) --compiler=$(CCNAME) --f90flags="$(CFLAGS)" \
		only: kernel_gaussian3d fft ifft conv_hist3d conv_prod conv_density : \
		$(FFTWFLAGS) \
		> /dev/null

readtree.so: readtree.f90
	$(F2PY) -m $(basename $@) -c $^ \
		--fcompiler=$(FCNAME) --compiler=$(CCNAME) --f90flags="$(CFLAGS)" \
		> /dev/null

extrema.pyf: extrema.f90
	$(F2PY) -h $@  --overwrite-signature -m $(basename $@) $^ skip: extrema_compute_ext : \
		> /dev/null

extrema.so: types.o convolution.o extrema.pyf extrema/extrema_types.f90 extrema/extrema_mod.f90 extrema.f90
	cd extrema && make objects
	$(F2PY) -m $(basename $@) -c $^ skip: extrema_compute_ext_ext : \
		-I./extrema/include -I./modules \
		$(MKLFLAGS) \
		--fcompiler=$(FCNAME) --compiler=$(CCNAME) \
		--f90flags="$(CFLAGS)" --link-lapack_opt \
		-L./ -lconvolution \
		> /dev/null

%.o: %.f90
	$(FCC) $(FFTWFLAGS) $(CFLAGS) -I../dist/mod -c $^ -o $@


clean:
	cd extrema && make cleanall
	rm -f *.so *.pyc *.pyf *.o
