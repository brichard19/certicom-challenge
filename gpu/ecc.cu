#include "ecmath.cuh"
#include "ecc.cuh"

#define get_global_id() (blockDim.x * blockIdx.x + threadIdx.x)

#define get_global_size() (blockDim.x * gridDim.x)

__device__ void print_big_int(uint131_t& x)
{
  printf("%.16lX%.16lX%.16lX", x.v[2], x.v[1], x.v[0]);
}

// Set all public keys to point-at-infinity
extern "C" __global__ void clear_public_keys(uint131_t* x, uint131_t* y, int count)
{
  int idx = get_global_id();
  int dim = get_global_size();
  
  for(int i = idx; i < count; i += dim) {
    set_point_at_infinity(x[i]);
  }
}

extern "C" __global__ void reset_counters(uint64_t* start_pos, uint64_t value, int count)
{
  int idx = get_global_id();
  int dim = get_global_size();
  
  for(int i = idx; i < count; i += dim) {
    start_pos[i] = value;
  }
}

__device__ int get_bit(uint131_t x, int bit)
{
  return (int)(x.v[bit /64] >> (bit % 64)) & 1;
}


// If the private key bit for P is 1, then add Q to P
extern "C" __global__ void do_step(uint131_t* global_px, uint131_t* global_py,
                                   uint131_t* global_rx, uint131_t* global_ry,
                                   uint131_t* mbuf, int count,
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
  uint131_t inverse = _one;


  // Perform Qx - Px and then multiply them together
  for(; i < count; i+=dim) {
    uint131_t px = global_px[i];

    if(result != NULL && (px.v[0] & dpmask) == 0) {
      // Record distinguished point
      int idx = atomicAdd(result_count, 1);

      DPResult r;

      r.a = priv_key_a[i];
      r.x = px;
      r.y = global_py[i];
      r.length = counter - start_pos[i];

      result[idx] = r;

      // Grab a new point from the staging buffer
      idx = atomicSub(staging_count, 1) - 1;
    
      uint131_t new_x = staging[idx].x;
      uint131_t new_y = staging[idx].y;
      priv_key_a[i] = staging[idx].a;

      start_pos[i] = counter;
      px = new_x;
      global_px[i] = new_x;
      global_py[i] = new_y;
    }

    // TODO: Proper mask
    int idx = px.v[0] & rmask;

    uint131_t rx = global_rx[idx];

    // Point addition, rx - px
    uint131_t t = sub(rx, px);

    inverse = mul(inverse, t);
    mbuf[i] = inverse;
  }

  // Perform inversion
  inverse = inv(inverse);

  // Start at last element (undo final loop counter add)
  i -= dim;

  // Complete addition

  for(; i >= gid; i-=dim) {
    uint131_t px = global_px[i];
    uint131_t py = global_py[i];

    int idx = px.v[0] & rmask;

    uint131_t rx = global_rx[idx];
    uint131_t ry = global_ry[idx]; 

    uint131_t s;

    if(i > gid) {

      // Get the 2nd-last element (product of all factors up to that number)
      // e.g. abcd
      uint131_t m = mbuf[i - dim];

      // Multiply to cancel out all factors except the last one
      // e.g. abcd * (abcde)^-1 = e^-1
      s = mul(inverse, m);

      // Cancel out from the inverse
      // e.g. abcde * e^-1 = abcd
      uint131_t diff = sub(rx, px);
      
      inverse = mul(inverse, diff);

    } else {
      s = inverse;
    }

    // Perform point addition
    uint131_t rise = sub(ry, py); 
    s = mul(s, rise);
    uint131_t s2 = square(s);

    uint131_t tmp1 = sub(s2, px);
    uint131_t x = sub(tmp1, rx);

    uint131_t tmp2 = sub(px, x);
    uint131_t tmp3 = mul(s, tmp2);
    uint131_t y = sub(tmp3, py);

    global_px[i] = x;
    global_py[i] = y;
  }
}




// If the private key bit for P is 1, then add Q to P
__device__ void batch_multiply_step(uint131_t* global_px, uint131_t* global_py, uint131_t* private_keys, uint131_t* global_qx, uint131_t* global_qy, uint131_t* mbuf, int priv_key_bit, int count)
{

  int gid = get_global_id();
  int dim = get_global_size();

  const uint131_t qx = global_qx[priv_key_bit];
  const uint131_t qy = global_qy[priv_key_bit];

  int i = gid;
  uint131_t inverse = _one;

  // Perform Qx - Px and then multiply them together
  for(; i < count; i+=dim) {

    uint131_t px = global_px[i];
    uint131_t t;

    int bit = get_bit(private_keys[i], priv_key_bit);

    if(!bit || is_infinity(px) || equal(px, qx)) {
      
      // Nothing to do, just use 1
      t = _one;
    } else {
      // Point addition, qx - px
      t = sub(qx, px);
    }

    inverse = mul(inverse, t);
    mbuf[i] = inverse;
  }

  // Perform inversion
  inverse = inv(inverse);

  // Start at last element (undo final loop counter add)
  i -= dim;

  // Complete addition

  for(; i >= gid; i-=dim) {
    uint131_t px = global_px[i];
    uint131_t py = global_py[i];

    int bit = get_bit(private_keys[i], priv_key_bit);

    if(!bit) {
      continue;
    } else if(is_infinity(px)) {
      global_px[i] = qx;
      global_py[i] = qy;

      continue;
    } else if(equal(px, qx)) {
      global_px[i] = global_qx[priv_key_bit + 1];
      global_py[i] = global_qy[priv_key_bit + 1];
      continue;
    }
    
    uint131_t s;
    if(i > gid) {

      // Get the 2nd-last element (product of all factors up to that number)
      // e.g. abcd
      uint131_t m = mbuf[i - dim];

      // Multiply to cancel out all factors except the last one
      // e.g. abcd * (abcde)^-1 = e^-1
      s = mul(inverse, m);

      // Cancel out from the inverse
      // e.g. abcde * e^-1 = abcd
      uint131_t diff = sub(qx, px);
      
      inverse = mul(inverse, diff);

    } else {
      s = inverse;
    }

    // Perform point addition
    uint131_t rise = sub(qy, py);
    s = mul(s, rise);
    uint131_t s2 = square(s);

    uint131_t tmp1 = sub(s2, px);
    uint131_t x = sub(tmp1, qx);

    uint131_t tmp2 = sub(px, x);
    uint131_t tmp3 = mul(s, tmp2);
    uint131_t y = sub(tmp3, py);

    global_px[i] = x;
    global_py[i] = y;
  }
}


extern "C" __global__ void batch_multiply(uint131_t* global_px, uint131_t* global_py, uint131_t* private_keys, uint131_t* mbuf, uint131_t* gx, uint131_t* gy, int priv_key_bit, int count)
{
    batch_multiply_step(global_px, global_py, private_keys, gx, gy, mbuf, priv_key_bit, count);
}

extern "C" __global__ void sanity_check(uint131_t* global_px, uint131_t* global_py, int count, int* errors)
{
  int gid = get_global_id();
  int dim = get_global_size();

  for(int i = gid; i < count; i += dim) {
    uint131_t x = global_px[i];
    uint131_t y = global_py[i];

    if(point_exists(x, y) == false) {
      atomicAdd(errors, 1);
    }
  }
}