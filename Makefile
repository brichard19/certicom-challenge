
TARGET_PLATFORMS ?= nvidia amd

# Default to all available processors; honor explicit jobs and inherited jobservers.
ifeq ($(MAKELEVEL),0)
# GNU make can omit environment job flags from MAKEFLAGS while parsing.
ifeq ($(filter -j% --jobs% --jobserver%,$(MAKEFLAGS) $(shell printf '%s' "$$MAKEFLAGS")),)
MAKEFLAGS += -j$(shell nproc)
endif
endif

VALID_PLATFORMS = amd nvidia
INVALID_PLATFORMS := $(filter-out $(VALID_PLATFORMS),$(TARGET_PLATFORMS))

# Fail early if invalids exist
ifneq ($(strip $(INVALID_PLATFORMS)),)
$(error Invalid TARGET_PLATFORMS: $(INVALID_PLATFORMS). Valid options are: $(VALID_PLATFORMS))
endif

CXX=g++

DEFINES=
CFLAGS=
INCLUDE=

# Debugging
ifeq ($(DEBUG),1)
CFLAGS+=-DDEBUG -g
else
CFLAGS+=-O2
endif

# MPI
ifeq ($(BUILD_MPI),1)
CXX=mpic++
CFLAGS+=-DBUILD_MPI
INCLUDE+=-I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi
endif

# Check for CUDA_HOME
ifneq ($(origin CUDA_HOME), environment)
CUDA_HOME=/usr/local/cuda
endif

# Check for ROCM_HOME
ifneq ($(origin ROCM_HOME), environment)
ROCM_HOME=/opt/rocm
endif

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
HIPCC=$(ROCM_HOME)/bin/hipcc

CFLAGS+=-std=c++20

# Directories
CUR_DIR:=$(CURDIR)
BUILD_DIR=$(CUR_DIR)
LIB_DIR=$(CUR_DIR)/lib
BIN_DIR=$(CUR_DIR)/bin
ROCM_LIB=$(ROCM_HOME)/lib
OBJDIR=$(CUR_DIR)/obj

INCLUDE+=-I$(CUR_DIR) -I$(CUR_DIR)/src/include -I$(CUR_DIR)/src/gpu -I$(CUR_DIR)/third_party/fmt/include

ifeq ($(BUILD_MPI),1)
MPI_LIBS+=-lmpi -lmpi_cxx
MPI_LINKER+=-L/usr/lib/x86_64-linux-gnu/openmpi/lib
endif



LINKER_RHO=-lfmt

# Keep backend-independent objects shared by tools, benchmarks, and tests.
MATH_SOURCES := ecc montgomery uint131 util
COMMON_OBJECTS := $(addprefix $(OBJDIR)/common/,$(addsuffix .o,ec_rho $(MATH_SOURCES)))
MATH_OBJECTS := $(addprefix $(OBJDIR)/common/,$(addsuffix .o,$(MATH_SOURCES)))
CPU_OBJECTS := $(addprefix $(OBJDIR)/common/,CPUPointFinder.o CPUPointFinderF2N.o)
FMT_OBJECTS := $(addprefix $(OBJDIR)/fmt/,format.o os.o)
FMT_LIBRARY := $(LIB_DIR)/libfmt.a

CPU_FLAGS = -DBUILD_CPU
AMD_FLAGS = $(CXX_CFLAGS_AMD) -DBUILD_GPU $(ROCM_INCLUDE)
NVIDIA_FLAGS = $(CXX_CFLAGS_NVIDIA) -DBUILD_GPU $(ROCM_INCLUDE) $(NVIDIA_INCLUDE)
HOST_FLAGS = $(CPPFLAGS) $(CFLAGS) $(INCLUDE) -Isrc
DEPFLAGS = -MMD -MP -MF $(@:.o=.d) -MT $@

TARGETS = tests benchmark_cpu rho_cpu
ifneq ($(filter nvidia,$(TARGET_PLATFORMS)),)
TARGETS += benchmark_nvidia rho_nvidia
endif
ifneq ($(filter amd,$(TARGET_PLATFORMS)),)
TARGETS += benchmark_amd rho_amd
endif

.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.PHONY: all third_party tests clean FORCE gpu_amd gpu_nvidia \
        benchmark_cpu benchmark_amd benchmark_nvidia rho_cpu rho_amd rho_nvidia rho_db rho_solve

all: $(TARGETS) rho_db rho_solve

# Preserve the existing command names; real output files control rebuilds.
third_party: $(FMT_LIBRARY)
gpu_amd: $(OBJDIR)/ecc_amd.co
gpu_nvidia: $(OBJDIR)/ecc_nvidia.co
benchmark_cpu: benchmark-cpu
benchmark_amd: benchmark-amd
benchmark_nvidia: benchmark-nvidia
rho_cpu: rho-cpu
rho_amd: rho-amd
rho_nvidia: rho-nvidia
rho_db: rho-db
rho_solve: rho-solve
tests: tests/math_tests tests/cpu_point_finder_tests

# Changing compilers or flags must invalidate cached objects. Keep the timestamp
# unchanged on ordinary builds so the dependency graph remains incremental.
CONFIG := $(OBJDIR)/build-config
CONFIG_VARS := CXX CPPFLAGS CFLAGS INCLUDE CPU_FLAGS AMD_FLAGS NVIDIA_FLAGS HIPCC \
               HIPCC_CFLAGS_AMD HIPCC_CFLAGS_NVIDIA CUDA_HOME ROCM_HOME \
               LDFLAGS LDLIBS LINKER_AMD LIBS_AMD LINKER_NVIDIA LIBS_NVIDIA \
               LINKER_RHO MPI_LINKER MPI_LIBS AR
shell_quote = '$(subst ','"'"',$(1))'
$(CONFIG): FORCE
	@mkdir -p $(@D)
	@printf '%s\n' $(foreach var,$(CONFIG_VARS),$(call shell_quote,$(var)=$($(var)))) > $@.tmp
	@cmp -s $@.tmp $@ && rm -f $@.tmp || mv -f $@.tmp $@

$(OBJDIR)/common/%.o: src/%.cpp $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CXX) $(HOST_FLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/cpu/%.o: src/%.cpp $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CXX) $(HOST_FLAGS) $(CPU_FLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/amd/%.o: src/%.cpp $(CONFIG) Makefile
	@mkdir -p $(@D)
	HIP_PLATFORM=amd $(CXX) $(HOST_FLAGS) $(AMD_FLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/nvidia/%.o: src/%.cpp $(CONFIG) Makefile
	@mkdir -p $(@D)
	HIP_PLATFORM=nvidia $(CXX) $(HOST_FLAGS) $(NVIDIA_FLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/tests/%.o: tests/%.cpp $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CXX) $(HOST_FLAGS) $(DEPFLAGS) -c $< -o $@

$(OBJDIR)/fmt/%.o: third_party/fmt/src/%.cc $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CXX) $(HOST_FLAGS) $(DEPFLAGS) -c $< -o $@

$(FMT_LIBRARY): $(FMT_OBJECTS)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

$(OBJDIR)/ecc_amd.co: src/gpu/ecc.cu $(CONFIG) Makefile
	@mkdir -p $(@D)
	HIP_PLATFORM=amd $(HIPCC) $(HIPCC_CFLAGS_AMD) -D__HIP_PLATFORM_AMD__ -Isrc -Isrc/gpu -Isrc/include -MMD -MP -MF $(@:.co=.d) -MT $@ -c $< -o $@

$(OBJDIR)/ecc_nvidia.co: src/gpu/ecc.cu $(CONFIG) Makefile
	@mkdir -p $(@D)
	HIP_PLATFORM=nvidia $(HIPCC) $(HIPCC_CFLAGS_NVIDIA) -D__HIP_PLATFORM_NVIDIA__ -Isrc -Isrc/gpu $(NVIDIA_INCLUDE) -Isrc/include -MMD -MP -MF $(@:.co=.d) -MT $@ -c $< -o $@

benchmark-cpu: $(OBJDIR)/cpu/benchmark.o $(CPU_OBJECTS) $(COMMON_OBJECTS) $(FMT_LIBRARY)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

rho-cpu: $(OBJDIR)/cpu/rho-main.o $(CPU_OBJECTS) $(COMMON_OBJECTS) $(FMT_LIBRARY)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(MPI_LINKER) $(MPI_LIBS) $(LDLIBS)

benchmark-amd: $(OBJDIR)/amd/benchmark.o $(OBJDIR)/amd/GPUPointFinder.o $(COMMON_OBJECTS) $(OBJDIR)/ecc_amd.co $(FMT_LIBRARY)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LINKER_AMD) $(LIBS_AMD) $(LDLIBS)

rho-amd: $(OBJDIR)/amd/rho-main.o $(OBJDIR)/amd/GPUPointFinder.o $(COMMON_OBJECTS) $(OBJDIR)/ecc_amd.co $(FMT_LIBRARY)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LINKER_AMD) $(LIBS_AMD) $(MPI_LINKER) $(MPI_LIBS) $(LDLIBS)

benchmark-nvidia: $(OBJDIR)/nvidia/benchmark.o $(OBJDIR)/nvidia/GPUPointFinder.o $(COMMON_OBJECTS) $(OBJDIR)/ecc_nvidia.co $(FMT_LIBRARY)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LINKER_NVIDIA) $(LIBS_NVIDIA) $(LDLIBS)

rho-nvidia: $(OBJDIR)/nvidia/rho-main.o $(OBJDIR)/nvidia/GPUPointFinder.o $(COMMON_OBJECTS) $(OBJDIR)/ecc_nvidia.co $(FMT_LIBRARY)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LINKER_NVIDIA) $(LIBS_NVIDIA) $(MPI_LINKER) $(MPI_LIBS) $(LDLIBS)

rho-db: $(OBJDIR)/common/rho-db.o $(COMMON_OBJECTS) $(FMT_LIBRARY)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

rho-solve: $(OBJDIR)/common/rho-solve.o $(COMMON_OBJECTS) $(FMT_LIBRARY)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

tests/math_tests: $(OBJDIR)/tests/math_tests.o $(MATH_OBJECTS)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

tests/cpu_point_finder_tests: $(OBJDIR)/tests/cpu_point_finder_tests.o $(CPU_OBJECTS) $(COMMON_OBJECTS)
	$(CXX) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

-include $(wildcard $(OBJDIR)/*/*.d $(OBJDIR)/*.d)

clean:
	$(RM) -r $(OBJDIR) $(LIB_DIR)
	$(RM) src/*.o third_party/fmt/*.o
	$(RM) rho-amd benchmark-amd rho-nvidia benchmark-nvidia rho-cpu benchmark-cpu
	$(RM) tests/math_tests tests/cpu_point_finder_tests rho-db rho-solve
