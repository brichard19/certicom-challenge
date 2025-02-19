#ifndef _HIP_HELPER_H
#define _HIP_HELPER_H

#include <hip/hip_runtime.h>
#include <sstream>

#define HIP_CALL(condition)\
{\
  hipError_t err = condition;\
  if(err != hipSuccess) {\
    std::stringstream ss;\
    ss << "GPU error " << err << " " << hipGetErrorString(err) << " " << __LINE__ << std::endl;\
    throw std::runtime_error(ss.str());\
  }\
}\

#define HIP_IGNORE(func) \
(void)func\

#endif