CXX=g++

DEFINES=
CFLAGS=
INCLUDE=

ifeq ($(CURVE_P131),1)
BIN_PREFIX:=p131-
DEFINES+=-DCURVE_P131
else ifeq ($(CURVE_P79),1)
BIN_PREFIX:=p79-
DEFINES+=-DCURVE_P79
else
# Default to P131
BIN_PREFIX:=p131-
DEFINES+=-DCURVE_P131
endif

# Debugging
ifeq ($(DEBUG),1)
DEFINES+=-DDEBUG -g
else
DEFINES+=-O2
endif

# MPI
ifeq ($(USE_MPI),1)
CXX=mpic++
CFLAGS+=-DUSE_MPI
INCLUDE+=-I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi
endif

TARGET_PLATFORMS=amd nvidia

CUDA_HOME=/usr/local/cuda
ROCM_HOME=/opt/rocm

# Nvidia target
export NVIDIA_INCLUDE=-I${CUDA_HOME}/include
export HIPCC_CFLAGS_NVIDIA=-Wno-deprecated-declarations
export CXX_CFLAGS_NVIDIA=-D__HIP_PLATFORM_NVIDIA__ -Wno-deprecated-declarations -Wno-return-local-addr
export LIBS_NVIDIA=-lcuda -lcudart
export LINKER_NVIDIA=-L${CUDA_HOME}/lib64 -L/mnt/c/Windows/System32/lxss/lib

# AMD target
export ROCM_INCLUDE=-I$(ROCM_HOME)/include
export HIPCC_CFLAGS_AMD=-Wno-deprecated-declarations -fPIE
export CXX_CFLAGS_AMD=-D__HIP_PLATFORM_AMD__
export LIBS_AMD=-lamdhip64
export LINKER_AMD=-L${ROCM_HOME}/lib

CFLAGS+=${DEFINES} -std=c++20

# Directories
CUR_DIR=$(shell pwd)
BUILD_DIR=$(CUR_DIR)
LIB_DIR=$(CUR_DIR)/lib
BIN_DIR=$(CUR_DIR)/bin
ROCM_LIB=$(ROCM_HOME)/lib

INCLUDE+=-I$(CUR_DIR) -I$(CUR_DIR)/include -I$(CUR_DIR)/gpu -I$(CUR_DIR)/third_party/json11 -I$(CUR_DIR)/third_party/fmt/include

ifeq ($(USE_MPI),1)
MPI_LIBS+=-lmpi -lmpi_cxx
MPI_LINKER+=-L/usr/lib/x86_64-linux-gnu/openmpi/lib
endif

export BUILD_DIR
export BIN_PREFIX
export LIB_DIR
export BIN_DIR
export INCLUDE
export CFLAGS

all:	rho benchmark tests tools

tools:	common_lib third_party
	make -C tools

rho:	main.cpp gpu_lib third_party common_lib
	mkdir -p ${BIN_DIR}
ifneq ($(filter nvidia, $(TARGET_PLATFORMS)),)
	HIP_PLATFORM=nvidia ${CXX} ${CFLAGS} ${CXX_CFLAGS_NVIDIA} main.cpp http_client.cpp ${INCLUDE} ${ROCM_INCLUDE} ${NVIDIA_INCLUDE} -o ${BIN_DIR}/${BIN_PREFIX}rho-nvidia -L${LIB_DIR}  ${MPI_LINKER} ${LINKER_NVIDIA} -L${ROCM_LIB} ${MPI_LIBS} -lgpupf-nvidia ${LIBS_NVIDIA} -ljson11 -lfmt -lcurl -lcommon
endif

ifneq ($(filter amd, $(TARGET_PLATFORMS)),)
	HIP_PLATFORM=amd ${CXX} ${CFLAGS} ${CXX_CFLAGS_AMD} main.cpp http_client.cpp ${INCLUDE} ${ROCM_INCLUDE} -o ${BIN_DIR}/${BIN_PREFIX}rho-amd -L${LIB_DIR} ${LINKER_AMD} -L${ROCM_LIB} ${MPI_LINKER} -lgpupf-amd ${LIBS_AMD} ${MPI_LIBS} -ljson11 -lfmt -lcurl -lcommon
endif

benchmark:	benchmark.cpp gpu_lib third_party
	mkdir -p ${BIN_DIR}
ifneq ($(filter nvidia, $(TARGET_PLATFORMS)),)
	HIP_PLATFORM=nvidia ${CXX} ${CFLAGS} ${CXX_CFLAGS_NVIDIA} benchmark.cpp ${INCLUDE} ${ROCM_INCLUDE} ${NVIDIA_INCLUDE} -o ${BIN_DIR}/${BIN_PREFIX}bench-nvidia -L${LIB_DIR} ${LINKER_NVIDIA} -L${ROCM_LIB} -lgpupf-nvidia ${LIBS_NVIDIA} -lfmt -lcommon
endif

ifneq ($(filter amd, $(TARGET_PLATFORMS)),)
	HIP_PLATFORM=amd ${CXX} ${CFLAGS} ${CXX_CFLAGS_AMD} benchmark.cpp ${INCLUDE} ${ROCM_INCLUDE} -o ${BIN_DIR}/${BIN_PREFIX}bench-amd -L${LIB_DIR} ${LINKER_AMD} -L${ROCM_LIB} -lgpupf-amd ${LIBS_AMD} -lfmt -lcommon
endif

tests:	common_lib
	make -C tests

gpu_lib: common_lib third_party
ifneq ($(filter nvidia, $(TARGET_PLATFORMS)),)
	HIP_PLATFORM=nvidia make -C gpu
endif

ifneq ($(filter amd, $(TARGET_PLATFORMS)),)
	HIP_PLATFORM=amd make -C gpu
endif

common_lib:	third_party
	make -C common

.PHONY: third_party
third_party:
	make -C third_party/fmt
	make -C third_party/json11

clean:
	make -C gpu clean
	make -C common clean
	make -C tests clean
	make -C third_party/fmt clean
	make -C third_party/json11 clean
	rm -rfv ${LIB_DIR} ${BIN_DIR} *.o
