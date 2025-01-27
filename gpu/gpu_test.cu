#include <hip/hip_runtime.h>
#include <random>
#include <iostream>
#include "ecmath.cuh"
#include "hip_helper.h"
#include <chrono>

#define TEST_SIZE 65536

double get_time()
{
  uint64_t ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

  return (double)ms / 1000.0;
}

__device__ void print_uint131(const uint131_t& x)
{
  int gid = blockDim.x * blockIdx.x + threadIdx.x;

  printf("%d %.16lx %.16lx %.16lx\n", gid, x.v[2], x.v[1], x.v[0]);
}


__constant__ uint131_t _r = {{0xd189e497ae0c5c29, 0x0cd212c781aea937, 0xe}};

__device__ uint131_t barrett(uint131_t a, uint131_t b)
{
  uint64_t ab[5];
  uint64_t abr[8] = {0};


  mul160x(a.v, b.v, ab); 

  // multiply by r
  for(int i = 0; i < 5; i++) {
    uint64_t high = 0;
    for(int j = 0; j < 3; j++) {
      uint128_t prod = abr[i + j] + (uint128_t)ab[i] * _r.v[j] + high;
      high = (uint64_t)(prod >> 64);
      abr[i + j] = (uint64_t)prod;
    }
    abr[i + 3] = high;
  }
  
  // printf("%.16lx %.16lx\n", abr[7], abr[6]);
  // printf("%.16lx %.16lx %.16lx\n", abr[5], abr[4], abr[3]);
  // printf("%.16lx %.16lx %.16lx\n", abr[2], abr[1], abr[0]);

  // Divide by 2^262
  abr[0] = abr[4] >> 6 | abr[5] << 58;
  abr[1] = abr[5] >> 6 | abr[6] << 58;
  abr[2] = abr[6] >> 6 | abr[7] << 58;
  abr[2] &= ((uint64_t)1 << 61) - 1;

  // printf("%.16lx %.16lx %.16lx\n", abr[2], abr[1], abr[0]);

  // Multiply by p
  uint64_t x[5];
  mul160x(abr, _p.v, x);
  //mont_mulModR(abr, _p.v, x);

  // Subtract from ab
  uint64_t diff = ab[0] - x[0];
  x[0] = diff;
  int borrow = diff > ab[0] ? 1 : 0;

  diff = ab[1] - x[1] - borrow;
  x[1] = diff;
  borrow = diff > ab[1] ? 1 : 0;

  diff = ab[2] - x[2] - borrow;
  x[2] = diff;
  borrow = diff > ab[2] ? 1 : 0;

  // Mod P
  bool gte = true;
  for(int i = 2; i >= 0; i--) {
    if(x[i] > _p.v[i]) {
      break;
    } else if(x[i] < _p.v[i]) {
      gte = false;
      break;
    }
  }

  if(gte) {
    int borrow = 0;

    uint64_t diff = x[0] - _p.v[0];
    borrow = diff > x[0] ? 1 : 0;
    x[0] = diff;
    
    diff = x[1] - _p.v[1] - borrow;
    borrow = diff > x[1] ? 1 : 0;
    x[1] = diff;
    
    diff = x[2] - _p.v[2] - borrow;
    x[2] = diff;
  }

  uint131_t z;
  for(int i = 0; i < 3; i++) {
    z.v[i] = x[i];
  }

  return z;
}

__global__ void barrett_perf_test(uint131_t* a, uint131_t* c, int n)
{
  int gid = blockIdx.x * blockDim.x + threadIdx.x;
  int grid_size = blockDim.x * gridDim.x;

  for(int i = gid; i < n; i += grid_size) {
    uint131_t x = barrett(a[i], a[i]);
    for(int j = 0; j < 128; j++) {
      x = barrett(a[i], x);
    }

    c[i] = x;
  }
}

__global__ void montgomery_perf_test(uint131_t* a, uint131_t* c, int n)
{
  int gid = blockIdx.x * blockDim.x + threadIdx.x;
  int grid_size = blockDim.x * gridDim.x;

  for(int i = gid; i < n; i += grid_size) {
    uint131_t x = mul(a[i], a[i]);
    for(int j = 0; j < 128; j++) {
      x = mul(a[i], x);
    }

    c[i] = x;
  }
}

template<int CURVE> __device__ void kernel_sub_test_impl(const uint131_t* x, const uint131_t* y, int n)
{
  int gid = blockDim.x * blockIdx.x + threadIdx.x;
  int grid_size = blockDim.x * gridDim.x;

  for(int i = gid; i < n; i += grid_size) {
    uint131_t s = sub<CURVE>(x[i], y[i]);
    uint131_t q = add(s, y[i]);

    uint131_t s2 = sub<CURVE>(y[i], x[i]);
    uint131_t q2 = add(s2, x[i]);

    if(equal(x[i], q) == false || equal(y[i], q2) == false) {
      print_uint131(x[i]);
      print_uint131(y[i]);
      print_uint131(s);
      print_uint131(q);
      print_uint131(s2);
      print_uint131(q2);
    }
    
    assert(equal(x[i], q));
    assert(equal(y[i], q2));
  }
}

extern "C" __global__ void kernel_sub_test(const uint131_t* x, const uint131_t* y, int n)
{
#if defined(CURVE_P131)
  kernel_sub_test_impl<131>(x, y, n);
#elif defined(CURVE_P79)
  kernel_sub_test_impl<79>(x, y, n);
#else
#error "Curve not defined"
#endif
}

template<int CURVE> __device__ void kernel_mul_test_impl(const uint131_t* x, const uint131_t* y, uint131_t* z, int n)
{
  int gid = blockDim.x * blockIdx.x + threadIdx.x;
  int grid_size = blockDim.x * gridDim.x;

  // Test algebraic property (x - y) * z = xz - yz
  for(int i = gid; i < n; i += grid_size) {
    uint131_t ls = mul(sub<CURVE>(x[i], y[i]), z[i]);

    uint131_t rs1 = mul(x[i], z[i]);
    uint131_t rs2 = mul(y[i], z[i]);
    uint131_t rs = sub<CURVE>(rs1, rs2);

    //assert(equal(ls, rs));
    if(equal(ls, rs) == false) {
      printf("Error:\n");
      print_uint131(ls);

      print_uint131(rs1);
      print_uint131(rs2);
      print_uint131(rs);
      assert(false);
    }
  }
}

extern "C" __global__ void kernel_mul_test(const uint131_t* x, const uint131_t* y, uint131_t* z, int n)
{
#if defined(CURVE_P131)
  kernel_mul_test_impl<131>(x, y, z, n);
#elif defined(CURVE_P79)
  kernel_mul_test_impl<79>(x, y, z, n);
#else
#error "Curve not defined"
#endif
}

template<int CURVE> __device__ void kernel_barrett_mul_test_impl(const uint131_t* x, const uint131_t* y, uint131_t* z, int n)
{
  int gid = blockDim.x * blockIdx.x + threadIdx.x;
  int grid_size = blockDim.x * gridDim.x;

  // Test algebraic property (x - y) * z = xz - yz
  for(int i = gid; i < n; i += grid_size) {
    uint131_t ls = barrett(sub<CURVE>(x[i], y[i]), z[i]);

    uint131_t rs1 = barrett(x[i], z[i]);
    uint131_t rs2 = barrett(y[i], z[i]);
    uint131_t rs = sub<CURVE>(rs1, rs2);

    //assert(equal(ls, rs));
    if(equal(ls, rs) == false) {
      printf("Error:\n");
      print_uint131(ls);

      print_uint131(rs1);
      print_uint131(rs2);
      print_uint131(rs);
      assert(false);
    }
  }
}

extern "C" __global__ void kernel_barrett_mul_test(const uint131_t* x, const uint131_t* y, uint131_t* z, int n)
{
#if defined(CURVE_P131)
  kernel_barrett_mul_test_impl<131>(x, y, z, n);
#elif defined(CURVE_P79)
  kernel_barrett_mul_test_impl<79>(x, y, z, n);
#else
#error "Curve not defined"
#endif

}

// Test algebraic property (a - b)(a + b) = a^2 - b^2
template<int CURVE> __device__ void kernel_square_test_impl(const uint131_t* x, const uint131_t* y, int n)
{
  int gid = blockDim.x * blockIdx.x + threadIdx.x;
  int grid_size = blockDim.x * gridDim.x;

  // Test algebraic property (x - y) * z = xz - yz
  for(int i = gid; i < n; i += grid_size) {
    uint131_t ls = mul(sub<CURVE>(x[i], y[i]), add(x[i], y[i]));
    uint131_t rs = sub<CURVE>(square(x[i]), square(y[i]));
    assert(equal(ls, rs));
  }
}

extern "C" __global__ void kernel_square_test(const uint131_t* x, const uint131_t* y, int n)
{
#if defined(CURVE_P131)
  kernel_square_test_impl<131>(x, y, n);
#elif defined(CURVE_P79)
  kernel_square_test_impl<79>(x, y, n);
#else
#error "Curve not defined"
#endif
}

template<int CURVE> __device__ void kernel_inv_test_impl(const uint131_t* x, int n)
{
  int gid = blockDim.x * blockIdx.x + threadIdx.x;
  int grid_size = blockDim.x * gridDim.x;

  // Test algebraic property (x + y) * z = xz + yz
  for(int i = gid; i < n; i += grid_size) {
    uint131_t k = x[i]; 
    uint131_t inverse = inv<CURVE>(k);
    uint131_t prod = mul(k, inverse);

    assert(equal(prod, _one));
  }
}

extern "C" __global__ void kernel_inv_test(const uint131_t* x, int n)
{
#if defined(CURVE_P131)
  kernel_inv_test_impl<131>(x, n);
#elif defined(CURVE_P79)
  kernel_inv_test_impl<79>(x, n);
#else
#error "Curve not defined"
#endif
}

class IntRNG {

private:
  std::random_device _rd;
  std::mt19937 _gen;
  std::uniform_int_distribution<uint64_t> _d;

public:
  IntRNG()
  {
    _gen = std::mt19937(_rd());
  }

  uint64_t next()
  {
    return _d(_gen);
  }
};

// Internal RNG
IntRNG _rng;


uint131_t gen_key()
{
  uint131_t k;

  for(int i = 0; i < 3; i++) {
    k.v[i] = _rng.next();
  }
  
#if defined(CURVE_P131)
    k.v[2] %= 5;
    k.v[1] %= 0x8e1d43f293469e33;
#elif defined(CURVE_P79)
    k.v[2] = 0;
    k.v[1] %= (0x62ce + 1);
    k.v[0] %= 0x5177412aca899cf5;
#else
#error "Curve is undefined"
#endif

  return k;
}


template<class... T> hipError_t hipLaunchKernel(void* kernel, dim3 gridDim, dim3 blockDim, size_t sharedMem, T... args)
{
  std::vector<void*> ptr = {&args...};

  return hipLaunchKernel(kernel, gridDim, blockDim, ptr.data(), sharedMem, 0);
}

bool mul_test()
{
  int n = TEST_SIZE;
  uint131_t* x_dev = nullptr;
  uint131_t* y_dev = nullptr;
  uint131_t* z_dev = nullptr;

  std::vector<uint131_t> z_host(n);


  HIP_CALL(hipMallocManaged(&x_dev, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&y_dev, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&z_dev, sizeof(uint131_t) * n));

  for(int i = 0; i < n; i++) {
    x_dev[i] = gen_key();
    y_dev[i] = gen_key();
    z_dev[i] = gen_key();
  }

  std::cout << "MUL test" << std::endl;
  HIP_CALL(hipLaunchKernel((void*)kernel_mul_test, dim3(8, 1, 1), dim3(32, 1, 1), 0, x_dev, y_dev, z_dev, n));
  HIP_CALL(hipDeviceSynchronize());

  HIP_CALL(hipFree(x_dev));
  HIP_CALL(hipFree(y_dev));
  HIP_CALL(hipFree(z_dev));

  return true;
}


bool mul_barrett_test()
{
  int n = TEST_SIZE;
  uint131_t* x_dev = nullptr;
  uint131_t* y_dev = nullptr;
  uint131_t* z_dev = nullptr;

  std::vector<uint131_t> z_host(n);


  HIP_CALL(hipMallocManaged(&x_dev, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&y_dev, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&z_dev, sizeof(uint131_t) * n));

  for(int i = 0; i < n; i++) {
    x_dev[i] = gen_key();
    y_dev[i] = gen_key();
    z_dev[i] = gen_key();
  }

  std::cout << "Barrett MUL test" << std::endl;
  HIP_CALL(hipLaunchKernel((void*)kernel_barrett_mul_test, dim3(8, 1, 1), dim3(32, 1, 1), 0, x_dev, y_dev, z_dev, n));
  HIP_CALL(hipDeviceSynchronize());

  HIP_CALL(hipFree(x_dev));
  HIP_CALL(hipFree(y_dev));
  HIP_CALL(hipFree(z_dev));

  return true;
}

bool inv_test()
{
  int n = TEST_SIZE;
  uint131_t* x_dev = nullptr;

  HIP_CALL(hipMallocManaged(&x_dev, sizeof(uint131_t) * n));

  for(int i = 0; i < n; i++) {
    x_dev[i] = gen_key();
  }

  std::cout << "INV test" << std::endl;
  HIP_CALL(hipLaunchKernel((void*)kernel_inv_test, dim3(8, 1, 1), dim3(32, 1, 1), 0, x_dev, n));
  HIP_CALL(hipDeviceSynchronize());

  HIP_CALL(hipFree(x_dev));

  return true;
}


bool square_test()
{
  int n = TEST_SIZE;
  uint131_t* x_dev = nullptr;
  uint131_t* y_dev = nullptr;

  std::vector<uint131_t> z_host(n);


  HIP_CALL(hipMallocManaged(&x_dev, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&y_dev, sizeof(uint131_t) * n));

  for(int i = 0; i < n; i++) {
    x_dev[i] = gen_key();
    y_dev[i] = gen_key();
  }

  std::cout << "SQR test" << std::endl;
  HIP_CALL(hipLaunchKernel((void*)kernel_square_test, dim3(8, 1, 1), dim3(32, 1, 1), 0, x_dev, y_dev, n));
  HIP_CALL(hipDeviceSynchronize());

  HIP_CALL(hipFree(x_dev));
  HIP_CALL(hipFree(y_dev));

  return true;
}

bool sub_test()
{
  int n = TEST_SIZE;
  uint131_t* x_dev = nullptr;
  uint131_t* y_dev = nullptr;

  std::vector<uint131_t> z_host(n);


  HIP_CALL(hipMallocManaged(&x_dev, sizeof(uint131_t) * n));
  HIP_CALL(hipMallocManaged(&y_dev, sizeof(uint131_t) * n));

  for(int i = 0; i < n; i++) {
    x_dev[i] = gen_key();
    y_dev[i] = gen_key();
  }

  std::cout << "SUB test" << std::endl;
  HIP_CALL(hipLaunchKernel((void*)kernel_sub_test, dim3(8, 1, 1), dim3(32, 1, 1), 0, x_dev, y_dev, n));
  HIP_CALL(hipDeviceSynchronize());

  HIP_CALL(hipFree(x_dev));
  HIP_CALL(hipFree(y_dev));
  
  return true;
}

std::string get_gpu_name(int device_id)
{
  hipDeviceProp_t props;
 
  HIP_CALL(hipGetDeviceProperties(&props, device_id));

  return std::string(props.name);
}



bool mul_perf_test()
{
  int n = 64 * 1024 * 1024;
  uint131_t* dev_a = nullptr;
  uint131_t* dev_c = nullptr;

  HIP_CALL(hipMalloc(&dev_a, sizeof(uint131_t) * n));
  HIP_CALL(hipMalloc(&dev_c, sizeof(uint131_t) * n));

  hipDeviceProp_t props;
  HIP_CALL(hipGetDeviceProperties(&props, 0));

  int blocks = props.multiProcessorCount * 8;

  printf("Barrett:\n");
  double t0 = get_time();
  for(int i = 0; i < 3; i++) {
    HIP_CALL(hipLaunchKernel((void*)barrett_perf_test, dim3(blocks, 1, 1), dim3(32, 1, 1), 0, dev_a, dev_c, n));
  }
  HIP_CALL(hipDeviceSynchronize());
  double t1 = get_time();
  printf("%f seconds\n", t1 - t0);

  printf("Montgomery:\n");
  t0 = get_time();
  for(int i = 0; i < 3; i++) {
    HIP_CALL(hipLaunchKernel((void*)montgomery_perf_test, dim3(blocks, 1, 1), dim3(32, 1, 1), 0, dev_a, dev_c, n));
  }
  HIP_CALL(hipDeviceSynchronize());
  t1 = get_time();
  printf("%f seconds\n", t1 - t0);

  HIP_CALL(hipFree(dev_a));
  HIP_CALL(hipFree(dev_c));

  return true;
}

int main(int argc, char**argv)
{
  bool pass = true;

  std::string name = get_gpu_name(0);

  std::cout << name << std::endl;

  pass &= sub_test();
  pass &= mul_test();
  pass &= mul_barrett_test();
  pass &= square_test();
  pass &= inv_test();
  pass &= mul_perf_test();

  return 0;
}