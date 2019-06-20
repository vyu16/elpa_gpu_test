FC     = ifort
FFLAGS = -O3

CUDA_LIB = -L/usr/local/cuda-10.0/lib64 -lcublas -lcudart
ELPA_VER = elpa-2019.05.001.rc1
ELPA_INC = -I/home/wy29/opt/$(ELPA_VER)/build/include/$(ELPA_VER)/modules
ELPA_LIB = -L/home/wy29/opt/$(ELPA_VER)/build/lib -lelpa
MATH_LIB = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

SRC = test_serial_real.f90 test_serial_cmplx.f90
EXE = $(SRC:.f90=.x)

all: $(EXE)

%.x: %.f90
	$(FC) -o $@ $< $(ELPA_INC) $(ELPA_LIB) $(CUDA_LIB) $(MATH_LIB)

clean:
	rm -f *.x
