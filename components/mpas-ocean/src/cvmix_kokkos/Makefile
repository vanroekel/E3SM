KOKKOS_PATH = kokkos
EXE_NAME="test_code"

default: build

OS := $(shell uname -s)
#NVCC_RESULT := $(shell which nvcc 2> NULL)
#NVCC_TEST := $(notdir $(NVCC_RESULT))
#ifeq ($(NVCC_TEST),nvcc)
#KOKKOS_DEVICES="Cuda"
#CXX := ${KOKKOS_PATH}/bin/nvcc_wrapper
#EXE = ${EXE_NAME}.cuda
#KOKKOS_ARCH = "BDW,Pascal60"
#KOKKOS_CUDA_OPTIONS = "enable_lambda"
#CPPFLAGS := -DGPU -lineinfo
#else
#KOKKOS_DEVICES="Serial"
#CXX := g++
#EXE = ${EXE_NAME}.host
#KOKKOS_ARCH = "SNB"
#endif

#CPPFLAGS += -DLMD94_OBL
#DEBUG = -g
#ifeq "$(Debug)" "true"
#CPPFLAGS += -DDEBUG_OUTPUT
#CXXFLAGS = -g -Wall -c $(CPPFLAGS)
#else
#CXXFLAGS = -O3 -Wall -c $(CPPFLAGS)
#endif

#F90=gfortran
#F90FLAGS= -Wno-character-truncation

#LDFLAGS =
LINK = ${F90}
#LIB = 

EXTRA_INC = 

DEPFLAGS = -M

include $(KOKKOS_PATH)/Makefile.kokkos

.PHONY: print_vars

print_vars:
	echo $(KOKKOS_LDFLAGS)

OBJ=cvmix2_compute_mixing.o vertical_mixing.o readinput.o profiles.o grid.o cvmix_lightweight_driver.o

build: $(EXE)

$(EXE): $(OBJ) libcvmix2.a 
	$(CXX) $(KOKKOS_LDFLAGS) -o $(EXE) $(OBJ) -L. -lmixing -lgfortran -lkokkos -ldl

libs: libcvmix2.a

libcvmix2.a: cvmix2_compute_mixing.o $(KOKKOS_LINK_DEPENDS)
ifeq ($(OS),Linux)
	ar cr -o libcvmix2.a $(KOKKOS_LINK_DEPENDS) *.o 
endif
ifeq ($(OS),Darwin)
	libtool -static -o libcvmix2.a $(KOKKOS_LINK_DEPENDS) *.o 
endif
	ranlib libcvmix2.a


#	$(LINK) $(KOKKOS_LDFLAGS) $(OBJ)  $(LDFLAGS) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean
	rm -f *.o *.mod $(EXE) libcvmix2.a

.cpp.o: 
	$(CXX) -c $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) $<

%.o: %.f90
	$(F90) -c $(DEBUG) $(F90FLAGS) $<

