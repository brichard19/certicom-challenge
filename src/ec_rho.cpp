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

  for(int i = 0; i < 32; i++) {
    a.push_back(ecc::genkey(rng));
  }
  for(int i = 0; i < 32; i++) {
    b.push_back(ecc::genkey(rng));
  }

  for(int i = 0; i < 32; i++) {
    RWPoint rw;
    rw.a = a[i];
    rw.b = b[i];
    auto p1 = ecc::mul(a[i], ecc::g());
    auto p2 = ecc::mul(b[i], ecc::q());
    rw.p = ecc::add(p1, p2);

    assert(ecc::exists(rw.p));
    rw_vec.push_back(rw);
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
      {"ecp131", ecc::ecpoint_t({{0x991395bdb3af97ba, 0xc74c6e35add57bf0, 0x00000002}},
                                {{0xf1fd45ecd69da948, 0xf8292e66fdd75b41, 0x00000001}})},
      {"ecp89", ecc::ecpoint_t({{0xa76f5de5725addf6, 0x000000000019c961, 0x00000000}},
                               {{0x801887d44f65fbc6, 0x00000000003d6ba5, 0x00000000}})},
      {"ecp79", ecc::ecpoint_t({{0xe02a1c4e11d34a22, 0x0000000000000608, 0x00000000}},
                               {{0x50ee4ba255488040, 0x0000000000002f1d, 0x00000000}})},
  };

  if(!ecc::is_equal(expected[name], sum)) {
    printf("ERROR: CHECKSUM FAILED\n");
    assert(1);
  }

  return rw_vec;
}

EncodedDP encode_dp(const DistinguishedPoint& dp)
{
  EncodedDP encoded;

  // Remove distinguished bits by shifting right then
  // copy into array
  uint131_t x2 = mont::rshift(dp.p.x, dp.dp_bits);

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

  assert(curve == 79 || curve == 131 || curve == 89);

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
  p.x = mont::lshift(p.x, dpbits);

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

std::vector<DistinguishedPoint> decode_dps(const uint8_t* bytes, size_t size, bool verify)
{
  BinaryDecoder decoder(bytes, size);
  DPHeader header = decoder.decode<DPHeader>();

  assert(header.curve_id == 79 || header.curve_id == 131);

  std::vector<DistinguishedPoint> dps;

  int payload_size = size - sizeof(DPHeader);
  assert(payload_size % header.count == 0);

  int field_size = payload_size / header.count;

  EncodedDP encoded;
  for(int i = 0; i < header.count; i++) {

    decoder.decode(&encoded, sizeof(encoded));
    dps.push_back(decode_dp(encoded, header.dp_bits, verify));
  }

  return dps;
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