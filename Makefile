TARGET_PLATFORMS=nvidia amd
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
export NVIDIA_INCLUDE=-I${CUDA_HOME}/include
export HIPCC_CFLAGS_NVIDIA=-Wno-deprecated-declarations -Xptxas -v
export CXX_CFLAGS_NVIDIA=-D__HIP_PLATFORM_NVIDIA__ -Wno-deprecated-declarations -Wno-return-local-addr
export LIBS_NVIDIA=-lcuda -lcudart
LINKER_NVIDIA=-L${CUDA_HOME}/lib64

# For WSL systems
LINKER_NVIDIA+=-L/mnt/c/Windows/System32/lxss/lib
export LINKER_NVIDIA


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
OBJDIR=$(CUR_DIR)/obj

INCLUDE+=-I$(CUR_DIR) -I$(CUR_DIR)/include -I$(CUR_DIR)/gpu -I$(CUR_DIR)/third_party/json11 -I$(CUR_DIR)/third_party/fmt/include

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

CPP_RHO := main.cpp GPUPointFinder.cpp  ec_rho.cpp  ecc.cpp  http_client.cpp  montgomery.cpp  uint131.cpp  util.cpp
CPP_RHO := $(addprefix src/, $(CPP_RHO))
CPP_BENCH := benchmark.cpp GPUPointFinder.cpp  ec_rho.cpp  ecc.cpp montgomery.cpp  uint131.cpp  util.cpp
CPP_BENCH := $(addprefix src/, $(CPP_BENCH))


all:	benchmark_nvidia rho_nvidia benchmark_amd rho_amd

.PHONY: third_party
third_party:
	make -C third_party/fmt
	make -C third_party/json11

gpu_nvidia:
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=nvidia hipcc -c src/ecc.cu -o $(OBJDIR)/ecc_nvidia.o $(HIPCC_CFLAGS_NVIDIA) -D__HIP_PLATFORM_NVIDIA__ -Isrc -I/usr/local/cuda/include -Isrc/include

gpu_amd:
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=amd hipcc -c src/ecc.cu -o $(OBJDIR)/ecc_amd.o $(HIPCC_CFLAGS_AMD) -D__HIP_PLATFORM_AMD__ -Isrc -Isrc/include

benchmark_nvidia:	third_party gpu_nvidia
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=nvidia g++ $(CPP_BENCH) $(OBJDIR)/ecc_nvidia.o -o benchmark-nvidia $(CXX_CFLAGS_NVIDIA) -D__HIP_PLATFORM_NVIDIA__ -Isrc -Isrc/include -Isrc -L$(LIB_DIR) -L$(ROCM_LIB) $(ROCM_INCLUDE) $(NVIDIA_INCLUDE) $(INCLUDE) $(LINKER_NVIDIA) $(LIBS_NVIDIA) -lfmt

rho_nvidia:	third_party gpu_nvidia
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=nvidia g++ $(CPP_RHO) $(OBJDIR)/ecc_nvidia.o -o rho-nvidia $(CXX_CFLAGS_NVIDIA) -D__HIP_PLATFORM_NVIDIA__ -Isrc -Isrc/include -Isrc -L$(LIB_DIR) -L$(ROCM_LIB) $(ROCM_INCLUDE) $(NVIDIA_INCLUDE) $(INCLUDE) $(LINKER_NVIDIA) $(LIBS_NVIDIA) -lfmt -ljson11 -lcurl

benchmark_amd:	third_party	gpu_amd
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=amd g++ $(CPP_BENCH) $(OBJDIR)/ecc_amd.o -o benchmark-amd -D__HIP_PLATFORM_AMD__ -Isrc -Isrc/include -Isrc -L$(LIB_DIR) -L$(ROCM_LIB) $(ROCM_INCLUDE) $(INCLUDE) $(LIBS_AMD) $(LINKER_AMD) -lfmt

rho_amd:	third_party gpu_amd
	mkdir -p $(OBJDIR)
	HIP_PLATFORM=amd g++ $(CPP_RHO) $(OBJDIR)/ecc_amd.o -o rho-amd -D__HIP_PLATFORM_AMD__ -Isrc -Isrc/include -Isrc -L$(LIB_DIR) -L$(ROCM_LIB) $(ROCM_INCLUDE) $(INCLUDE) $(LINKER_AMD) $(LIBS_AMD) -lfmt -ljson11 -lcurl

clean:
	rm -rf obj