#include <cassert>
#include <stdexcept>

#include "ec_rho.h"

#include "rw_p79.h"
#include "rw_p131.h"
#include "montgomery.h"

#include "binary_encoder.h"

std::vector<RWPoint> get_rw_points()
{
  std::vector<RWPoint> rw_vec;

  for(int i = 0; i < 32; i++) {
    RWPoint rw;

    if(ecc::curve_name() == "ecp131") {
      rw.a = make_uint131(_p131_a_str[i]);
      rw.b = make_uint131(_p131_b_str[i]);
      rw.p = ecc::ecpoint_t(make_uint131(_p131_x_str[i]), make_uint131(_p131_y_str[i]));
    } else if(ecc::curve_name() == "ecp79") {
      rw.a = make_uint131(_p79_a_str[i]);
      rw.b = make_uint131(_p79_b_str[i]);
      rw.p = ecc::ecpoint_t(make_uint131(_p79_x_str[i]), make_uint131(_p79_y_str[i]));
    } else {
      throw std::runtime_error("No curve selected");
    }

    // Validate
    ecc::ecpoint_t p = ecc::add(ecc::mul(rw.a, ecc::g()), ecc::mul(rw.b, ecc::q()));

    assert(ecc::is_equal(p, rw.p));

    rw_vec.push_back(rw);
  }

  return rw_vec;
}



// Encodes a DistinguishedPoint into a string of bytes
std::vector<uint8_t> encode_dp(const DistinguishedPoint& dp, int dbits)
{
  std::vector<uint8_t> buf(128);
  uint8_t* ptr = buf.data();

  // Remove distinguished bits by shifting right then
  // copy into array
  uint131_t x2 = mont::rshift(dp.p.x, dbits);
  
  int len = ((131 - dbits) + 7) / 8;
  memcpy(ptr, &x2, len);
  ptr += len;

  // sign bit: 1 byte
  uint8_t sign = is_odd(dp.p.y) ? 1 : 0;
  *ptr = sign;
  ptr++;

  // exponent
  len = (131 + 7) / 8;
  memcpy(ptr, &dp.a, len);
  ptr += len;

  // Walk length: 40 bits (5 bytes)
  memcpy(ptr, &dp.length, 5);
  ptr += 5;

  // Resize the array
  buf.resize(ptr - buf.data());

  return buf;
}

std::vector<uint8_t> encode_dps(const std::vector<DistinguishedPoint>& dps, int dbits, int curve)
{

  assert(curve == 79 || curve == 131);

  BinaryEncoder encoder;

  // Encode header
  DPHeader header;
  header.version = 1;
  header.count = dps.size();
  header.dbits = dbits;
  header.curve = curve;

  encoder.encode(header);

  // Encode points
  for(auto dp : dps) {
    auto buf = encode_dp(dp, dbits);
    encoder.encode(buf.data(), buf.size());
  }

  // Convert bytes to vector
  std::vector<uint8_t> vec(encoder.get_size());
  memcpy(vec.data(), encoder.get_ptr(), encoder.get_size());

  return vec;
}

DistinguishedPoint decode_dp(const uint8_t* bytes, int dbits)
{
  ecc::ecpoint_t p;
  memset(&p, 0, sizeof(p));

  const uint8_t* ptr = bytes;

  const int a_bytes = (131 + 7) / 8;
  const int x_bytes = ((131 - dbits) + 7) / 8;

  // extract x
  memcpy(&p.x, ptr, x_bytes);
  p.x = mont::lshift(p.x, dbits);
  ptr += x_bytes;

  // sign
  uint8_t sign = *ptr;
  ptr++;
 
  // Calculate y component
  p.y = ecc::calc_y(p.x, sign);

  uint131_t a;
  memset(&a, 0, sizeof(a));
  memcpy(&a, ptr, a_bytes);
  ptr += a_bytes;

  uint64_t length = 0;
  memcpy(&length, ptr, 5);

  assert(ecc::exists(p));

  return DistinguishedPoint(a, p, length);
}

std::vector<DistinguishedPoint> decode_dps(const uint8_t* bytes, size_t size)
{
  BinaryDecoder decoder(bytes, size);
  DPHeader header = decoder.decode<DPHeader>();
  
  assert(header.curve == 79 || header.curve == 131);

  std::vector<DistinguishedPoint> dps;

  int payload_size = size - sizeof(DPHeader);
  assert(payload_size % header.count == 0);

  int field_size = payload_size / header.count;
  std::vector<uint8_t> buf(field_size);

  for(int i = 0; i < header.count; i++) {
    decoder.decode(buf.data(), buf.size());

    dps.push_back(decode_dp(buf.data(), header.dbits));
  }

  return dps;
}