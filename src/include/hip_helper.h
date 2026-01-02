#ifndef _HIP_HELPER_H
#define _HIP_HELPER_H

#include <hip/hip_runtime.h>
#include <sstream>
#include <string>
#include <vector>

struct HIPDeviceID {
  int id;
  std::string uuid;

  HIPDeviceID(int id, std::string uuid) : id(id), uuid(uuid) {}
};

using HIPDeviceMap = std::vector<HIPDeviceID>;

#define HIP_CALL(condition)                                                                        \
  {                                                                                                \
    hipError_t err = condition;                                                                    \
    if(err != hipSuccess) {                                                                        \
      std::stringstream ss;                                                                        \
      ss << "GPU error " << err << " " << hipGetErrorString(err) << " " << __FILE__ << ":"         \
         << __LINE__ << std::endl;                                                                 \
      throw std::runtime_error(ss.str());                                                          \
    }                                                                                              \
  }

#define HIP_IGNORE(func) (void)func

template <class... T>
hipError_t hipLaunchKernel(void* kernel, dim3 gridDim, dim3 blockDim, size_t sharedMem, T... args)
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

#ifdef __HIP_PLATFORM_AMD__
  // TODO: Find a way to check for WGP mode. Seems hipcc defaults to WGP mode for RDNA.
  std::string arch = std::string(props.gcnArchName);
  if(arch.find("gfx10") != std::string::npos || arch.find("gfx11") != std::string::npos ||
     arch.find("gfx12") != std::string::npos) {
    return props.multiProcessorCount * 2;
  } else {
    return props.multiProcessorCount;
  }
#else
  return props.multiProcessorCount;
#endif
}

inline int get_warp_size(int device_id)
{
  hipDeviceProp_t props;

  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return props.warpSize;
}

inline int get_compute_mode(int device_id)
{
  hipDeviceProp_t props;

  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return props.computeMode;
}

inline HIPDeviceMap get_device_map()
{
  int count;

  HIP_CALL(hipGetDeviceCount(&count));

  HIPDeviceMap device_map;

  for(int d = 0; d < count; d++) {
    hipDeviceProp_t props;

    HIP_CALL(hipGetDeviceProperties(&props, d));

    std::string uuid(props.uuid.bytes, 16);

    device_map.push_back(HIPDeviceID(d, uuid));
  }

  return device_map;
}

#endif