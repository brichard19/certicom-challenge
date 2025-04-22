
ifeq ($(HIP_PLATFORM), nvidia)
LIB_NAME=libgpupf-nvidia.a
HIPCC_CFLAGS=${HIPCC_CFLAGS_NVIDIA} --save-temps -g
GPU_LINKER=${LINKER_NVIDIA}
GPU_LIBS=${LIBS_NVIDIA}
ROCM_INCLUDE:=${ROCM_INCLUDE} ${NVIDIA_INCLUDE}
CFLAGS:=${CFLAGS}
CXX_CFLAGS:=${CXX_CFLAGS_NVIDIA}
else
LIB_NAME=libgpupf-amd.a
HIPCC_CFLAGS=${HIPCC_CFLAGS_AMD}
CFLAGS:=${CFLAGS}
GPU_LINKER=${LINKER_AMD}
GPU_LIBS=${LIBS_AMD}
CXX_CFLAGS:=${CXX_CFLAGS_AMD}
endif


cuda_lib:		GPUPointFinder.cpp ecc.cu
	g++ ${CFLAGS} ${CXX_CFLAGS} -c GPUPointFinder.cpp ${ROCM_INCLUDE} ${INCLUDE} -o GPUPointFinder.o
	hipcc -c ecc.cu -o ecc_cu.o ${HIPCC_CFLAGS} ${CFLAGS} ${ROCM_INCLUDE} ${INCLUDE}
	mkdir -p ${LIB_DIR}
	ar rcs ${LIB_DIR}/${LIB_NAME} GPUPointFinder.o ecc_cu.o

ifeq ($(HIP_PLATFORM), nvidia)
	HIP_PLATFORM=nvidia hipcc gpu_test.cu -o gpu_test-nvidia ${CFLAGS} ${ROCM_INCLUDE} ${INCLUDE} -L${LIB_DIR} ${GPU_LINKER} -L${ROCM_LIB} ${GPU_LIBS}
endif

ifeq ($(HIP_PLATFORM), amd)
	HIP_PLATFORM=amd hipcc gpu_test.cu -o gpu_test-amd ${CFLAGS} ${ROCM_INCLUDE} ${INCLUDE} -L${LIB_DIR} ${GPU_LINKER} -L${ROCM_LIB} ${GPU_LIBS}
endif

clean:
	rm -rfv *.o gpu_test-amd gpu_test-nvidia