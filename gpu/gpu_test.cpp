#include <assert.h>
#include <hip/hip_runtime.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdint.h>

#include "ec_rho.h"
#include "hip_helper.h"
#include "mont.h"
#include "util.h"

namespace {


std::string get_gpu_name(int device_id)
{
  hipDeviceProp_t props;
 
  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return std::string(props.name);
}

template<class... T> hipError_t hipLaunchKernel(void* kernel, dim3 gridDim, dim3 blockDim, size_t sharedMem, T... args)
{
  std::vector<void*> ptr = {&args...};

  return hipLaunchKernel(kernel, gridDim, blockDim, ptr.data(), sharedMem, 0);
}


}

extern "C" __global__ void kernel_mul_test(const uint131_t* x, const uint131_t* y, uint131_t* z, int n);


bool mul_test()
{
  int n = 1024;
  uint131_t* x_dev = nullptr;
  uint131_t* y_dev = nullptr;
  uint131_t* z_dev = nullptr;

  std::vector<uint131_t> z_host(n);


  HIP_CALL(hipMallocManaged(&x_dev, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&y_dev, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&z_dev, sizeof(uint131_t) * n));

  std::cout << "Generating keys.." << std::endl;
  for(int i = 0; i < n; i++) {
    x_dev[i] = ecc::genkey();
    y_dev[i] = ecc::genkey();
  }

  std::cout << "Launching kernel" << std::endl;
  HIP_CALL(hipLaunchKernel((void*)kernel_mul_test, dim3(1, 1, 1), dim3(32, 1, 1), 0, x_dev, y_dev, z_dev, n));
  HIP_CALL(hipDeviceSynchronize());

  std::cout << "Done" << std::endl;

  for(int i = 0; i < n; i++) {
    z_host[i] = mont::mul(x_dev[i], y_dev[i]);
  }

  for(int i = 0; i < n; i++) {
    if(z_dev[i] != z_host[i]) {
      std::cout << "FAIL" << std::endl;

      std::cout << to_str(x_dev[i]) << std::endl;
      std::cout << to_str(y_dev[i]) << std::endl;
      std::cout << to_str(z_dev[i]) << std::endl;
      std::cout << to_str(z_host[i]) << std::endl;

      return false;
    }
  }

  HIP_CALL(hipFree(x_dev));
  HIP_CALL(hipFree(y_dev));
  HIP_CALL(hipFree(z_dev));

  return true;
}

int main(int argc, char**argv)
{
#if defined(CURVE_P131)
  ecc::set_curve("ecp131");
#elif defined(CURVE_P79)
  ecc::set_curve("ecp79");
#else
#error "Curve is undefined"
#endif

  bool pass = true;

  std::string name = get_gpu_name(0);

  std::cout << name << std::endl;

  pass &= mul_test();

  return 0;
}