#ifndef _HIP_HELPER_H
#define _HIP_HELPER_H

#include <hip/hip_runtime.h>
#include <sstream>
#include <vector>
#include <string>

#define HIP_CALL(condition)\
{\
  hipError_t err = condition;\
  if(err != hipSuccess) {\
    std::stringstream ss;\
    ss << "GPU error " << err << " " << hipGetErrorString(err) << " " << __FILE__ << ":" << __LINE__ << std::endl;\
    throw std::runtime_error(ss.str());\
  }\
}\

#define HIP_IGNORE(func) \
(void)func\

template<class... T> hipError_t hipLaunchKernel(void* kernel, dim3 gridDim, dim3 blockDim, size_t sharedMem, T... args)
{
  std::vector<void*> ptr = {&args...};

  return hipLaunchKernel(kernel, gridDim, blockDim, ptr.data(), sharedMem, 0);
}

inline std::string get_gpu_name(int device_id)
{
  hipDeviceProp_t props;
 
  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return std::string(props.name);
}

inline int get_cu_count(int device_id)
{
  hipDeviceProp_t props;
 
  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return props.multiProcessorCount;
}

inline int get_warp_size(int device_id)
{
  hipDeviceProp_t props;
 
  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return props.warpSize;
}

#endif