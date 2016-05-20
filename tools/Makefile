
F2PY=f2py
FCC=ifort
CCNAME=intelem
FCNAME=intelem

MODS=-I../dist/mod -module modules
INCLUDEDIR=../dist-no-p/mod

MKLARCH = em64t
MKLPATH=$(MKLROOT)/lib/$(MKLARCH)

LAPACK  = -I$(MKLROOT)/include -mkl
LAPACKLIB = $(LAPACK) -lpthread
CFLAGS=-Ofast -fp-model strict -static -free -cpp -xHost $(LAPACKLIB) -traceback -mkl -lpthread -fPIC

MKLFLAGS=-L/softs/fftw/3.3.4-pureintel/lib -I./include -I/softs/intel/compilers_and_libraries_2016.3.210/linux/mkl/include
FFTWFLAGS=-lfftw3 -L/softs/fftw/3.3.4-pureintel/lib -I/softs/fftw/3.3.4-pureintel/include

all: _tools.so convolution.so readtree.so extrema.so

_tools.pyf: types.f90 compute.f90 misc.f90 io.f90
	$(F2PY) -h $@ --overwrite-signature -m $(basename $@) $^ skip: parse_params read_info_headers \
		read_brick_header read_brick_data read_particle read_mergertree_headers_1 \
		read_mergertree_headers_2 read_mergertree_parent_of read_1 read_1_dynamics \
		read_2_dummy read_2 read_region :  \
		> /dev/null

_tools.so: types.o misc.o compute.o io.o _tools.pyf
	$(F2PY) -c $(LIBS) $^ -I$(INCLUDEDIR) -m $(basename $@)  skip: parse_params read_info_headers\
		read_brick_header read_brick_data read_particle read_mergertree_headers_1\
		read_mergertree_headers_2 read_mergertree_parent_of read_1 read_1_dynamics\
		read_2_dummy read_2 read_region : \
		--compiler=$(CCNAME) --fcompiler=$(FCNAME) --f90flags="$(CFLAGS)" \
		> /dev/null

convolution.so: types.o convolution.o
	$(F2PY) -m $(basename $@) -c $^ \
		--fcompiler=$(FCNAME) --compiler=$(CCNAME) --f90flags="$(CFLAGS)" \
		only: kernel_gaussian3d fft ifft conv_hist3d conv_prod conv_density : \
		$(FFTWFLAGS)\
		> /dev/null

readtree.so: readtree.f90
	$(F2PY) -m $(basename $@) -c $^ \
		--fcompiler=$(FCNAME) --compiler=$(CCNAME) --f90flags="$(CFLAGS)" \
		> /dev/null

extrema.pyf: extrema.f90
	$(F2PY) -h $@ -m $(basename $@) $^ --overwrite-signature > /dev/null

%.o: %.f90
	$(FCC) -c $(CFLAGS) $(FFTWFLAGS) $(MKLFLAGS) $^ -o $@ $(MODS)

%: %.so

extrema: make_extrema extrema.so

make_extrema:
	cd extrema && make objects
	cp extrema/include/* modules/ -v
extrema.so: convolution.o extrema.o extrema.pyf
	$(F2PY) -m $(basename $@) -c $^\
		skip: extrema_compute_ext : \
		-I./extrema/include \
		$(MKLFLAGS) $(FFTWFLAGS) \
		--fcompiler=$(FCNAME) --compiler=$(CCNAME) \
		--f90flags="$(CFLAGS)" --link-lapack_opt

clean:
	cd extrema && make cleanall
	rm -f *.mod *.o *.so *.pyc *.pyf
