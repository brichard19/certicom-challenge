#include <cassert>
#include <map>
#include <stdexcept>

#include "binary_encoder.h"
#include "ec_rho.h"
#include "montgomery.h"

std::vector<RWPoint> get_rw_points()
{
  std::vector<RWPoint> rw_vec;

  std::string name = ecc::curve_name();

  std::vector<uint131_t> a;
  std::vector<uint131_t> b;
  DeterministicRNG rng(0x1234);

  int count = 0;

  while(count < 32) {

    uint131_t a = ecc::genkey(rng);
    uint131_t b = ecc::genkey(rng);

    auto p1 = ecc::mul(a, ecc::g());
    auto p2 = ecc::mul(b, ecc::q());
    auto p = ecc::add(p1, p2);

    // We only want p where the high 3 bits of x and y are zero so that x and y both fit in 128 bits
    if((p.x.w.v2 == 0) && (p.y.w.v2 == 0)) {
      RWPoint rw;
      rw.a = a;
      rw.b = b;
      rw.p = p;

      assert(ecc::exists(rw.p));
      rw_vec.push_back(rw);
      count++;
    }
  }

  // Verify
  ecc::ecpoint_t sum;
  for(int i = 0; i < 32; i++) {
    sum = ecc::add(sum, rw_vec[i].p);
  }

  assert(ecc::exists(sum));

#if defined(DEBUG)
  printf("Checksum:\n");
  printf("{{0x%.16lx, 0x%.16lx, 0x%.8x}}\n", sum.x.w.v0, sum.x.w.v1, sum.x.w.v2);
  printf("{{0x%.16lx, 0x%.16lx, 0x%.8x}}\n", sum.y.w.v0, sum.y.w.v1, sum.y.w.v2);
#endif

  std::map<std::string, ecc::ecpoint_t> expected = {
      {"ec2n131", ecc::ecpoint_t({{0xa81d4846a24d3869, 0x7aaa68638b9ecc4e, 0x00000000}},
                                 {{0x4281ad579379762d, 0x81d07eaffd9a8e89, 0x00000007}})},

      {"ec2n89", ecc::ecpoint_t({{0x7c638e4ed3b0eef7, 0x00000000017f17f2, 0x00000000}},
                                {{0x789575e57d671fb4, 0x00000000002a3fab, 0x00000000}})},

      {"ec2n79", ecc::ecpoint_t({{0xec000b39c2afd195, 0x0000000000007f1b, 0x00000000}},
                                {{0x3d13df6ced2b1989, 0x0000000000003041, 0x00000000}})},

      {"ecp131", ecc::ecpoint_t({{0x35d0286229b66c14, 0x749c0f06a7121cbc, 0x00000000}},
                                {{0x5f4c3aa56f0fb31d, 0x3286b1bd2be15ecf, 0x00000001}})},

      {"ecp89", ecc::ecpoint_t({{0x17d5a9d43f7ae3ef, 0x00000000001e1c5d, 0x00000000}},
                               {{0xeef7866b17c91b37, 0x0000000000e637ce, 0x00000000}})},

      {"ecp79", ecc::ecpoint_t({{0x352235cf969c42c8, 0x00000000000057d4, 0x00000000}},
                               {{0x8c539fa7274f69a0, 0x00000000000055ec, 0x00000000}})},
  };

  auto expected_sum = expected.find(name);
  if(expected_sum != expected.end() && !ecc::is_equal(expected_sum->second, sum)) {
    printf("ERROR: CHECKSUM FAILED\n");
    assert(false);
  }

  return rw_vec;
}

EncodedDP encode_dp(const DistinguishedPoint& dp)
{
  EncodedDP encoded;

  // Remove distinguished bits by shifting right then
  // copy into array
  uint131_t x2 = rshift(dp.p.x, dp.dp_bits);

  memcpy(encoded.tx, &x2, sizeof(encoded.tx));

  // sign bit: 1 byte
  encoded.data.sign = is_odd(dp.p.y) ? 1 : 0;

  memcpy(encoded.data.a, &dp.a, sizeof(encoded.data.a));

  memcpy(encoded.len, &dp.length, sizeof(encoded.len));

  encoded.checksum = (uint8_t)dp.p.y.v[0];

  return encoded;
}

std::vector<uint8_t> encode_dps(const std::vector<DistinguishedPoint>& dps, int curve, int dpbits)
{

  BinaryEncoder encoder(dps.size() * sizeof(EncodedDP));

  // Encode header
  DPHeader header;
  header.version = 1;
  header.count = dps.size();
  header.curve_id = curve;
  header.dp_bits = dpbits;

  encoder.encode(header);

  // Encode points
  for(auto dp : dps) {
    EncodedDP encoded = encode_dp(dp);
    encoder.encode(&encoded, sizeof(encoded));
  }

  // Convert bytes to vector
  std::vector<uint8_t> vec(encoder.get_size());
  memcpy(vec.data(), encoder.get_ptr(), encoder.get_size());

  return vec;
}

DistinguishedPoint decode_dp(const EncodedDP& dp, int dpbits, bool verify)
{
  ecc::ecpoint_t p;
  memset(&p, 0, sizeof(p));

  // extract x
  memcpy(&p.x, dp.tx, sizeof(dp.tx));
  p.x = lshift(p.x, dpbits);

  // sign
  uint8_t sign = dp.data.sign;

  // Calculate y component
  p.y = ecc::calc_y(p.x, sign);

  uint131_t a;
  memset(&a, 0, sizeof(a));
  memcpy(&a, dp.data.a, sizeof(dp.data.a));

  uint64_t length = 0;
  memcpy(&length, dp.len, sizeof(dp.len));

  assert(ecc::exists(p));

  if(verify) {
    assert(dp.checksum == (uint8_t)p.y.v[0]);
  }
  return DistinguishedPoint(a, p, dpbits, length);
}

bool verify_dp(const DistinguishedPoint& dp)
{
  auto r_points = get_rw_points();

  auto p = ecc::mul(dp.a, ecc::g());
  auto key_a = dp.a;
  auto key_b = make_uint131(0);

  uint32_t mask = (1 << dp.dp_bits) - 1;
  for(uint64_t i = 0; i < dp.length; i++) {
    int idx = p.x.v[0] & 0x1f;

    p = ecc::add(p, r_points[idx].p);
    key_a = ecc::add_priv_keys(key_a, r_points[idx].a);
    key_b = ecc::add_priv_keys(key_b, r_points[idx].b);
  }

  ecc::ecpoint_t p2 = ecc::add(ecc::mul(key_a, ecc::g()), ecc::mul(key_b, ecc::q()));
  if(!ecc::is_equal(p, p2)) {
    return false;
  }

  if(!ecc::is_equal(p, dp.p)) {
    return false;
  }

  return true;
}
