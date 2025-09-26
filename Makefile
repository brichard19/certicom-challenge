
TARGET_PLATFORMS ?= nvidia amd

VALID_PLATFORMS = amd nvidia
INVALID_PLATFORMS := $(filter-out $(VALID_PLATFORMS),$(TARGET_PLATFORMS))

# Fail early if invalids exist
ifneq ($(strip $(INVALID_TARGETS)),)
$(error Invalid TARGETS: $(INVALID_PLATFORMS). Valid options are: $(VALID_PLATFORMS))
endif

CXX=g++

DEFINES=
CFLAGS=
INCLUDE=

# Debugging
ifeq ($(DEBUG),1)
DEFINES+=-DDEBUG -g
else
DEFINES+=-O2
endif

# MPI
ifeq ($(BUILD_MPI),1)
CXX=mpic++
CFLAGS+=-DBUILD_MPI
INCLUDE+=-I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi
endif

CUDA_HOME=/usr/local/cuda
ROCM_HOME=/opt/rocm

# Nvidia target
NVIDIA_INCLUDE=-I${CUDA_HOME}/include
HIPCC_CFLAGS_NVIDIA=-Wno-deprecated-declarations -Xptxas -v
CXX_CFLAGS_NVIDIA=-D__HIP_PLATFORM_NVIDIA__ -Wno-deprecated-declarations -Wno-return-local-addr -std=c++20
LIBS_NVIDIA=-lcuda -lcudart
LINKER_NVIDIA=-L${CUDA_HOME}/lib64

IS_WSL := $(shell uname -r | grep -i microsoft)

# For WSL systems
ifneq ($(IS_WSL),)
LINKER_NVIDIA+=-L/mnt/c/Windows/System32/lxss/lib
endif


# AMD target
ROCM_INCLUDE=-I$(ROCM_HOME)/include
HIPCC_CFLAGS_AMD=-Wno-deprecated-declarations -fPIE
CXX_CFLAGS_AMD=-D__HIP_PLATFORM_AMD__ -std=c++20
LIBS_AMD=-lamdhip64
LINKER_AMD=-L${ROCM_HOME}/lib

#CFLAGS+=${DEFINES} -std=c++20

# Directories
CUR_DIR=$(shell pwd)
BUILD_DIR=$(CUR_DIR)
LIB_DIR=$(CUR_DIR)/lib
BIN_DIR=$(CUR_DIR)/bin
ROCM_LIB=$(ROCM_HOME)/lib
OBJDIR=$(CUR_DIR)/obj

INCLUDE+=-I$(CUR_DIR) -I$(CUR_DIR)/src/include -I$(CUR_DIR)/gpu -I$(CUR_DIR)/third_party/json11 -I$(CUR_DIR)/third_party/fmt/include

ifeq ($(BUILD_MPI),1)
MPI_LIBS+=-lmpi -lmpi_cxx
MPI_LINKER+=-L/usr/lib/x86_64-linux-gnu/openmpi/lib
endif

export BUILD_DIR
export BIN_PREFIX=certicom-
export LIB_DIR
export BIN_DIR
export INCLUDE
export CFLAGS

CPP_MATH_TESTS := ecc.cpp montgomery.cpp uint131.cpp util.cpp
CPP_MATH_TESTS := $(addprefix src/, $(CPP_MATH_TESTS))
CPP_RHO := main.cpp GPUPointFinder.cpp  ec_rho.cpp  ecc.cpp  http_client.cpp  montgomery.cpp  uint131.cpp  util.cpp
CPP_RHO := $(addprefix src/, $(CPP_RHO))
CPP_BENCH := benchmark.cpp GPUPointFinder.cpp  ec_rho.cpp  ecc.cpp montgomery.cpp  uint131.cpp  util.cpp
CPP_BENCH := $(addprefix src/, $(CPP_BENCH))

TARGETS = tests

# NVIDIA targets
ifeq ($(filter nvidia,$(TARGET_PLATFORMS)),nvidia)
TARGETS += benchmark_nvidia rho_nvidia
endif

ifeq ($(filter amd,$(TARGET_PLATFORMS)),amd)
TARGETS += benchmark_amd rho_amd
endif

all:	$(TARGETS)

.PHONY: third_party
third_party:
	make -C third_party/json11

gpu_nvidia:
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=nvidia hipcc -c src/ecc.cu -o $(OBJDIR)/ecc_nvidia.co $(HIPCC_CFLAGS_NVIDIA) -D__HIP_PLATFORM_NVIDIA__ -Isrc -I/usr/local/cuda/include -Isrc/include

gpu_amd:
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=amd hipcc -c src/ecc.cu -o $(OBJDIR)/ecc_amd.co $(HIPCC_CFLAGS_AMD) -D__HIP_PLATFORM_AMD__ -Isrc -Isrc/include

benchmark_nvidia:	third_party gpu_nvidia
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=nvidia $(CXX) $(CPP_BENCH) $(OBJDIR)/ecc_nvidia.co -o benchmark-nvidia $(CXX_CFLAGS_NVIDIA) -D__HIP_PLATFORM_NVIDIA__ -Isrc -Isrc/include -Isrc -L$(LIB_DIR) -L$(ROCM_LIB) $(ROCM_INCLUDE) $(NVIDIA_INCLUDE) $(INCLUDE) $(LINKER_NVIDIA) $(LIBS_NVIDIA)

rho_nvidia:	third_party gpu_nvidia
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=nvidia $(CXX) $(CPP_RHO) $(OBJDIR)/ecc_nvidia.co -o rho-nvidia $(CXX_CFLAGS_NVIDIA) -D__HIP_PLATFORM_NVIDIA__ -Isrc -L$(LIB_DIR) -L$(ROCM_LIB) $(ROCM_INCLUDE) $(NVIDIA_INCLUDE) $(INCLUDE) $(LINKER_NVIDIA) $(LIBS_NVIDIA) -ljson11 -lcurl

benchmark_amd:	third_party	gpu_amd
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=amd $(CXX) $(CPP_BENCH) $(OBJDIR)/ecc_amd.co -o benchmark-amd $(CXX_CFLAGS_AMD) -Isrc -Isrc -L$(LIB_DIR) -L$(ROCM_LIB) $(ROCM_INCLUDE) $(INCLUDE) $(LIBS_AMD) $(LINKER_AMD)

rho_amd:	third_party gpu_amd
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=amd $(CXX) $(CPP_RHO) $(OBJDIR)/ecc_amd.co -o rho-amd $(CXX_CFLAGS_AMD) -Isrc -Isrc/include -Isrc -L$(LIB_DIR) -L$(ROCM_LIB) $(ROCM_INCLUDE) $(INCLUDE) $(LINKER_AMD) $(LIBS_AMD) -ljson11 -lcurl

.PHONY: tests
tests:
	$(CXX) tests/math_tests.cpp $(CPP_MATH_TESTS) -o tests/math_tests $(INCLUDE)

clean:
	rm -v -rf src/*.o
	rm -v -rf obj
	rm -v -f rho-amd
	rm -v -f benchmark-amd
	rm -v -f rho-nvidia
	rm -v -f benchmark-nvidia