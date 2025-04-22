#include "ecmath.cuh"
#include "ecc.cuh"

#define get_global_id() (blockDim.x * blockIdx.x + threadIdx.x)

#define get_global_size() (blockDim.x * gridDim.x)

__device__ void print_big_int(uint131_t& x)
{
  printf("%.8X%.16lX%.16lX", x.w.v2, x.w.v1, x.w.v0);
}


__device__ int get_bit(uint131_t x, int bit)
{
  if(bit >= 128) {
    return (x.w.v2 >> (bit & 0x7f)) & 1;
  } else if(bit >= 64) {
    return (x.w.v1 >> (bit & 0x3f)) & 1;
  } else {
    return (x.w.v0 >> (bit & 0x1f)) & 1;
  }
}

__device__ uint131_t load(const uint131_vec_ptr_t p, int n, int i)
{
  uint131_t value;

  __int128* p128 = (__int128*)p;
  uint32_t* p32 = (uint32_t*)p;
  
  __int128 lo = p128[i];
  uint32_t hi = p32[n * 4 + i];

  value.w.v0 = (uint64_t)lo;
  value.w.v1 = (uint64_t)(lo >> 64);
  value.w.v2 = hi;

  return value;
}

__device__ void store(uint131_vec_ptr_t p, int n, int i, uint131_t x)
{
  __int128* p128 = (__int128*)p;
  uint32_t* p32 = (uint32_t*)p;

  __int128 lo = (__int128)x.w.v1 << 64 | x.w.v0;
  p128[i] = lo;

  p32[4 * n + i] = x.w.v2;
}

// If the private key bit for P is 1, then add Q to P
template<int CURVE> __device__ void do_step_impl(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py,
                                   uint131_vec_ptr_t global_rx, uint131_vec_ptr_t global_ry,
                                   uint131_vec_ptr_t mbuf, int count,
                                   DPResult* result, int* result_count,
                                   StagingPoint* staging, int* staging_count,
                                   uint131_t* priv_key_a,
                                   uint64_t counter,
                                   uint64_t* start_pos,
                                   uint64_t dpmask)
{
  const int rmask = 0x1f;
  const int gid = get_global_id();
  const int dim = get_global_size();

  int i = gid;

  uint131_t inverse;

  // Perform Qx - Px and then multiply them together
  for(; i < count; i+=dim) {
    uint131_t px = load(global_px, count, i);

    // TODO: Proper mask
    int idx = px.w.v0 & rmask;

    uint131_t rx = load(global_rx, 32, idx);

    // Point addition, rx - px
    uint131_t t = sub<CURVE>(rx, px);
    if(i > gid) {
      inverse = mul<CURVE>(inverse, t);
    } else {
      inverse = t;
    }
    store(mbuf, count, i, inverse);
  }

  // Perform inversion
  inverse = inv<CURVE>(inverse);

  // Start at last element (undo final loop counter add)
  i -= dim;

  // Complete addition

  for(; i >= gid; i-=dim) {
    uint131_t px = load(global_px, count, i);
    uint131_t py = load(global_py, count, i);

    int idx = px.w.v0 & rmask;

    uint131_t rx = load(global_rx, 32, idx);
    uint131_t ry = load(global_ry, 32, idx);

    uint131_t s;

    if(i > gid) {

      // Get the 2nd-last element (product of all factors up to that number)
      // e.g. abcd
      uint131_t m = load(mbuf, count, i - dim);

      // Multiply to cancel out all factors except the last one
      // e.g. abcd * (abcde)^-1 = e^-1
      s = mul<CURVE>(inverse, m);

      // Cancel out from the inverse
      // e.g. abcde * e^-1 = abcd
      uint131_t diff = sub<CURVE>(rx, px);
      
      inverse = mul<CURVE>(inverse, diff);

    } else {
      s = inverse;
    }

    // Perform point addition
    uint131_t rise = sub<CURVE>(ry, py);
    s = mul<CURVE>(s, rise);
    uint131_t s2 = square<CURVE>(s);

    uint131_t tmp1 = sub<CURVE>(s2, px);
    uint131_t x = sub<CURVE>(tmp1, rx);

    uint131_t tmp2 = sub<CURVE>(px, x);
    uint131_t tmp3 = mul<CURVE>(s, tmp2);
    uint131_t y = sub<CURVE>(tmp3, py);

    store(global_px, count, i, x);
    store(global_py, count, i, y);
  }
}




// If the private key bit for P is 1, then add Q to P
template<int CURVE> __device__ void batch_multiply_step(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py, uint131_t* private_keys, uint131_t* global_qx, uint131_t* global_qy, uint131_vec_ptr_t mbuf, int priv_key_bit, int count)
{
  int gid = get_global_id();
  int dim = get_global_size();

  const uint131_t qx = global_qx[priv_key_bit];
  const uint131_t qy = global_qy[priv_key_bit];

  int i = gid;

  uint131_t one;
  if(CURVE == 131) {
    one = _p131_one;
  } else if(CURVE == 79) {
    one = _p79_one;
  }

  uint131_t inverse = one;
  
  // Perform Qx - Px and then multiply them together
  for(; i < count; i+=dim) {

    uint131_t px = load(global_px, count, i);
    uint131_t t;

    int bit = get_bit(private_keys[i], priv_key_bit);

    if(!bit || is_infinity(px) || equal(px, qx)) {
      
      // Nothing to do, just use 1
      t = one;
    } else {
      // Point addition, qx - px
      t = sub<CURVE>(qx, px);
    }

    inverse = mul<CURVE>(inverse, t);
    store(mbuf, count, i, inverse);
  }

  // Perform inversion
  inverse = inv<CURVE>(inverse);

  // Start at last element (undo final loop counter add)
  i -= dim;

  // Complete addition

  for(; i >= gid; i-=dim) {
    uint131_t px = load(global_px, count, i);
    uint131_t py = load(global_py, count, i);

    int bit = get_bit(private_keys[i], priv_key_bit);

    if(!bit) {
      continue;
    } else if(is_infinity(px)) {
      store(global_px, count, i, qx);
      store(global_py, count, i, qy);

      continue;
    } else if(equal(px, qx)) {
      store(global_px, count, i, load(global_qx, count, priv_key_bit + 1));
      store(global_py, count, i, load(global_qy, count, priv_key_bit + 1));
      continue;
    }
    
    uint131_t s;
    if(i > gid) {

      // Get the 2nd-last element (product of all factors up to that number)
      // e.g. abcd
      uint131_t m = load(mbuf, count, i - dim);

      // Multiply to cancel out all factors except the last one
      // e.g. abcd * (abcde)^-1 = e^-1
      s = mul<CURVE>(inverse, m);

      // Cancel out from the inverse
      // e.g. abcde * e^-1 = abcd
      uint131_t diff = sub<CURVE>(qx, px);
      
      inverse = mul<CURVE>(inverse, diff);

    } else {
      s = inverse;
    }

    // Perform point addition
    uint131_t rise = sub<CURVE>(qy, py);
    s = mul<CURVE>(s, rise);
    uint131_t s2 = square<CURVE>(s);

    uint131_t tmp1 = sub<CURVE>(s2, px);
    uint131_t x = sub<CURVE>(tmp1, qx);

    uint131_t tmp2 = sub<CURVE>(px, x);
    uint131_t tmp3 = mul<CURVE>(s, tmp2);
    uint131_t y = sub<CURVE>(tmp3, py);

    store(global_px, count, i, x);
    store(global_py, count, i, y);
  }
}

template<int CURVE> __device__ void sanity_check_impl(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py, int count, int* errors)
{
  int gid = get_global_id();
  int dim = get_global_size();

  for(int i = gid; i < count; i += dim) {
    uint131_t x = load(global_px, count, i);
    uint131_t y = load(global_py, count, i);

    if(point_exists<CURVE>(x, y) == false) {
      atomicAdd(errors, 1);
    }
  }
}

__device__ void clear_public_keys_impl(uint131_vec_ptr_t x, uint131_vec_ptr_t y, int count)
{
  int idx = get_global_id();
  int dim = get_global_size();
 
  uint131_t point_zero;
  set_point_at_infinity(point_zero);
  for(int i = idx; i < count; i += dim) {
    store(x, count, i, point_zero);
  }
}

__device__ void reset_counters_impl(uint64_t* start_pos, uint64_t value, int count)
{
  int idx = get_global_id();
  int dim = get_global_size();
  
  for(int i = idx; i < count; i += dim) {
    start_pos[i] = value;
  }
}

// Set all public keys to point-at-infinity
extern "C" __global__ void clear_public_keys(uint131_vec_ptr_t x, uint131_vec_ptr_t y, int count)
{
  clear_public_keys_impl(x, y, count);
}

extern "C" __global__ void reset_counters(uint64_t* start_pos, uint64_t value, int count)
{
  reset_counters_impl(start_pos, value, count);
}



extern "C" __global__ void sanity_check_p131(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py, int count, int* errors)
{
  sanity_check_impl<131>(global_px, global_py, count, errors);
}

extern "C" __global__ void sanity_check_p79(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py, int count, int* errors)
{
  sanity_check_impl<79>(global_px, global_py, count, errors);
}


extern "C" __global__ void batch_multiply_p79(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py, uint131_t* private_keys, uint131_vec_ptr_t mbuf, uint131_t* gx, uint131_t* gy, int priv_key_bit, int count)
{
    batch_multiply_step<79>(global_px, global_py, private_keys, gx, gy, mbuf, priv_key_bit, count);
}


extern "C" __global__ void do_step_p79(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py,
                                   uint131_vec_ptr_t global_rx, uint131_vec_ptr_t global_ry,
                                   uint131_vec_ptr_t mbuf, int count,
                                   DPResult* result, int* result_count,
                                   StagingPoint* staging, int* staging_count,
                                   uint131_t* priv_key_a,
                                   uint64_t counter,
                                   uint64_t* start_pos,
                                   uint64_t dpmask)
{
  do_step_impl<79>(global_px, global_py, global_rx, global_ry, mbuf, count, result, result_count, staging, staging_count, priv_key_a, counter, start_pos, dpmask);
}



extern "C" __global__ void batch_multiply_p131(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py, uint131_t* private_keys, uint131_vec_ptr_t mbuf, uint131_t* gx, uint131_t* gy, int priv_key_bit, int count)
{
    batch_multiply_step<131>(global_px, global_py, private_keys, gx, gy, mbuf, priv_key_bit, count);
}


extern "C" __global__ void check_for_dp(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py,
                                   int count,
                                   DPResult* result, int* result_count,
                                   StagingPoint* staging, int* staging_count,
                                   uint131_t* priv_key_a,
                                   uint64_t counter,
                                   uint64_t* start_pos,
                                   uint64_t dpmask)
{
  const int gid = get_global_id();
  const int dim = get_global_size();

  int i = gid;

  for(; i < count; i+=dim) {
    uint131_t px = load(global_px, count, i);
    if(result != NULL && (px.w.v0 & dpmask) == 0) {
      // Record distinguished point
      int idx = atomicAdd(result_count, 1);

      DPResult r;

      r.a = priv_key_a[i];
      r.x = px;
      r.y = load(global_py, count, i);
      r.length = counter - start_pos[i];

      result[idx] = r;

      // Grab a new point from the staging buffer
      idx = atomicSub(staging_count, 1) - 1;
    
      uint131_t new_x = staging[idx].x;
      uint131_t new_y = staging[idx].y;
      priv_key_a[i] = staging[idx].a;

      start_pos[i] = counter;
      px = new_x;
      store(global_px, count, i, new_x);
      store(global_py, count, i, new_y);
    }
  }
}

extern "C" __global__ void do_step_p131(uint131_vec_ptr_t global_px, uint131_vec_ptr_t global_py,
                                   uint131_vec_ptr_t global_rx, uint131_vec_ptr_t global_ry,
                                   uint131_vec_ptr_t mbuf, int count,
                                   DPResult* result, int* result_count,
                                   StagingPoint* staging, int* staging_count,
                                   uint131_t* priv_key_a,
                                   uint64_t counter,
                                   uint64_t* start_pos,
                                   uint64_t dpmask)
{
  do_step_impl<131>(global_px, global_py, global_rx, global_ry, mbuf, count, result, result_count, staging, staging_count, priv_key_a, counter, start_pos, dpmask);
}