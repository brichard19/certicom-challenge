#include <cassert>
#include <stdexcept>

#include "ec_rho.h"
#include "montgomery.h"
#include "binary_encoder.h"

namespace {
// DO NOT MODIFY
std::string _p79_a_str[] = {
"c806168f60b90d7c758",
"1984517914b97ccc1364",
"2deed643606e2a22002c",
"af0f8b9174d05c07ff7",
"52248cab7d2291be77c5",
"1eb1936a6303310fac6",
"3306f509824ae973bf13",
"20f62dcd7a058e7db9be",
"120f8c2c38be74942c23",
"338ae6f3737a41124f74",
"187af490e26e243b3342",
"205a50c11193447299e8",
"5bc57ded78e9916896b9",
"251b8826f8d9148bbda0",
"5b90888d2b82e1721d4f",
"333a7f456008f1b48773",
"3e8782c66d7308c3e084",
"54acb0faaec3c08319ad",
"359843717ee3e17ebf0b",
"4af891f0c1e7efbe66b8",
"2e895547a16c1087c764",
"41b522d8cda0ce5b924e",
"430d751c31f1009d6faf",
"27c9f8c726d3c43f2138",
"5cc55299d57f2a654270",
"49ec897ac98154b79d6c",
"651c04fac912c07c709",
"4487e2b97379c31598a6",
"48a3bf746dc70b4e662a",
"254577444d1878c51498",
"3e0fdfc769e22a80ed7c",
"3d06c06c898a002c96ea",
};

std::string _p79_b_str[] = {
"23573b8c0594b154a7fc",
"1211e94deb5d9f959750",
"5cadcdce5e6cac5dd968",
"600e7edc5eb8a2af3ce5",
"2fc433085787f95c8ccd",
"11e05f6e995f65e4a401",
"318e67b52daa9b357062",
"2b12f4fb6f1b65b77469",
"38b864a4a65da340606a",
"1e0b684473b6df8556b4",
"470b7bed297003c91919",
"226c78e7dd14d0cd3096",
"428a4e1fa65f3e35fbc2",
"1894e84e4e962cdfbb42",
"1859387671918e78902a",
"59399b32e3c2e7afda6a",
"41df7b3657c29c5e26e3",
"204ddc5ee834d5e81b3d",
"59c8a416e159b19f51d7",
"30353613c9028867a174",
"45a1ea9f070ae65d7f9a",
"a861f66444221e79dd4",
"3932db73ffdff4526882",
"ec3ac32591ae83e6558",
"69a23a39ef0de86e928",
"250b525ba49a1bc18810",
"108e6248774e4309b24f",
"2845de39f3531db60878",
"1f0ab0eb7640e52e89df",
"4d8d32889a552bf91898",
"500a3f277a391ef173ed",
"3162535bf0766f03758b",
};

std::string _p79_x_str[] = {
"5eea148f33b6c74b6b70",
"58556493bd797ea4dfb4",
"c6dc75dd16412af22b8",
"47cb84c1661cf0bca84a",
"caf3d75c132a1be0311",
"535fd64e5b5579a3cbf4",
"565c6f6a1af365f75450",
"4abbbc7bb0da83bca543",
"58a983871e3b80c19f5a",
"49c9927a486ada8d7c4b",
"193de36abec10ddda583",
"4b89b92ba3ca9a9791e4",
"1a80ecfebdc93c8d4b3f",
"5b866b72e65e341ac2bf",
"1ab0a878947080106ba1",
"4a02ac6e02b493031880",
"5c20f9ac7a24e81a89c8",
"5951c89372ff3be3f5c8",
"33436672a7ecb1425a0",
"122891921cc51f14cb04",
"23133a45a47502f4625d",
"1e94c0579d09134f4f2b",
"9cecee4e6b205cdbfcb",
"3200ae379fde01cac131",
"23e071a6f30730f7c617",
"5c79d22058adc0cbaadc",
"154ec7dcd68478621d2a",
"3a6addfa0f7963f58cda",
"353ed348a0c68c6876ba",
"20fea812c32439bddfb7",
"479709d014008b1b9f47",
"e80ea68fea1bfaa5929",
};

std::string _p79_y_str[] = {
"18f97ea75cc6c4627f60",
"3264807b895a41e61e6b",
"4c0a8d3fd894ee3c267b",
"5a92aef821158fe9f5e2",
"2036680eaeb1c40ae8b5",
"3a4cd519d9f884f37ded",
"39bbc7a02ea3b28bec48",
"5afafe764ab8be3437eb",
"4036895340f114492478",
"3a35e3ee2a06831c3698",
"cc39d120fbf61b985a9",
"3cc613b3daf3e94c422",
"1ede89597af4529b9461",
"11782b128a07c141a036",
"55f1ef546dd3a6a3b5f4",
"17ffa492b50bcc52ccec",
"4cffde78edc4555dd973",
"3fa04ef27c900e1d1aa7",
"144d72c52fb735494b0e",
"5bc6cc4a9041dde06540",
"5ade0cc26001b0ca24b1",
"b9e8da03c238fe7d49e",
"61910597b8defac7644b",
"22d7066b86e80f83d88f",
"27b6cc72c9b5980492c5",
"1b2dac6193b05dedf608",
"148a480b4b938a2bf903",
"4e5448d7d576f9643995",
"2c327ec80ebaec9f9c9d",
"599e540383e11577969b",
"d4ef05ff7548aae57d5",
"5704abe009a23a727e89",
};

// DO NOT MODIFY
std::string _p131_a_str[] = {
"395099fedd0c41ba333de03c47d276521",
"164c3f82247df3c898e06103f3b62a45",
"1a046381ad95f352c4e348480e6fab685",
"10e79625e05e2ae86b130837aff5298a1",
"2c590637af9b87acab9b099396ab02106",
"14385a85fa3f442b4670a8284e9feb95f",
"385947f69a81349195dabb361880a1fbd",
"38be60cc631145f34a775ae52249d7dd0",
"cef395abe17806220ad61e74f7f6e981",
"4458e5b071775057b5bf22b35d69cc9c2",
"2b5e33955610a64eefca33716200f757d",
"5a07c5a50caf4f9cbe95054a46fe611c",
"409232580e80df4257376a9b1f1e61401",
"18d05654a9b9e06751947cf0b5662ee00",
"197450a4cb7e9db88a7ff2478467aebc6",
"76b950843a657c011c52d65071ef4176",
"3d21d453336b4c9bc9efdf66ffa1d141e",
"22169cd12fb66a8d4c9188a8bfdce0e8c",
"217848576657f3b58923c1666080168e8",
"31f1e511059c4c2fbacf4538377e1bbff",
"52bfb40bffa618fde2769b6572afa360",
"2b362c9f09e2301f376d44673094cee79",
"29074606af56907c32f4fc2fe0b51f3cb",
"37f309703abda7981efadb45bb1e970b2",
"20c3cf99599a8fffb48c64ac1f8110f6c",
"253b477b49cff38611433746f6151e10",
"2bbbedc60f6ae4910fa69d2278b75d558",
"3f4efa61bf810b56461e60bec1bc5699c",
"21e33db10301d4af1e6e88f48eb3251e7",
"4e700aefe549d47f3e12c47e5960e06",
"32e35a05cb269806d1c491180a8a7c6b0",
"19b5439e55d8972bea7e72d6ee9b13305",
};

std::string _p131_b_str[] = {
"3d9753cb0ed05be729c2548147f0890a8",
"20996ddc7e51956bade3efab040da93a",
"20ebe0e330512617cdc5a3773bd30a449",
"424521a2c0f7ceabc8be22703d4203559",
"2637ce693fc2f5009e5f5ed369fb67aad",
"f4f24938f2a4a7344f1dcd976705e976",
"2e79d16049c34aa558bda966c538bba30",
"425f438316c27f3112c12b85709cc5e2e",
"3b66a77fbb1f4340d0feead8b0409cdf1",
"4319f7e5b2b37b92cd2de4665e5a73316",
"1cfad41299e44636b63c1e477f73e5e89",
"32629dedd1251524d8fc71be13b2bc5a5",
"3920c53196564498ba547be67bb380cce",
"2d7c3990c7df678702635a3b7210a1622",
"475ac4ae9ecc4c4ee154cd98aac57be57",
"418d551e838851a8c419c9a0feaf9171e",
"54a3c9c75e461ab203dc67beed2f2c63",
"a546120031ac0b2dc2788e6db04e0b13",
"543e9633c2436481e5cc09d31a763d84",
"86da6c96d39ce8233afdaeaa9c6323b1",
"337375ebdef0c5c28e1b04cc581e4dea0",
"28498c311494e6b8c9eeec26e9e00aed5",
"36f6917f9dd9529941cc9b1a559187ce0",
"113843c63d9feddd408bf42e4ee2ad5ed",
"47f94781e672550f6d88bf4a04a5a68d6",
"168bd295a90abd3fe5fd12536a8e316c",
"26c299f1dbd4ae3094f8474660518942b",
"4102b0f5b956eb2a98e929bca328cbc2a",
"18a490a8634a0937d95bd0f19f742486b",
"27e7b622826dc0825ba83b7d8d6c57494",
"12ad55fb6707800bed3577bf639c40251",
"3f8cd4892bd7fb92b2bfa037a303ee19f",
};

std::string _p131_x_str[] = {
"2f1a93815a28180514204b71d4b7a1d10",
"29a14adf24c6a9feb6a2b014fcda606aa",
"371f0eef6f3a3ab2f71a47e695eb5bbb7",
"43bd9147c72034f00d6ae1bdca47005a5",
"45bce51bcea32af479c0d3fd36983c723",
"3c2ebba543c6cf178beadd1ef90a42135",
"310cc64a34e3e9192d029a42bb22c5afc",
"16428e593979107754e4e1431c9fa0f3e",
"38d8ac1fe17d612f1721ee5c3b7f189a9",
"2c7757833dda7968ad9743edb0ae92cdc",
"232248e8ff12c16b57cd87ad0cc72d46a",
"465dd2bb82d89c3ab0098aff040470d6c",
"17485194c828b031c8cb19e8865be17f7",
"46dacc3cd6e98f9723cba3be325ccc63b",
"19dfb2e69a1b96a9fcc7b4e992f0b433f",
"17f6af59a820d4290946da459610968ff",
"3dd8396dac6f269a039bebb248348f2fd",
"1728eb8d72d60173ad8aeb2ef1d1bc7fa",
"36d4bdec4b07451f505f826bbd32c4bf6",
"2dfbc231ee0fdb40738ff46435f26833d",
"48177947fa7496a3892f6f6374cfc70e2",
"21edf989ba7e01fc03c508d73d46b9ca6",
"13ec7dd47a558f16ecaae87c3ea4e5cea",
"afa75fdc37467ca207afc8170cfcf8c5",
"a7cdc9338d18e24d7743bd4670b5fb1e",
"4073167245012c69691517586e4ea099c",
"438001ec2b9c87f1c52306a49499187fc",
"630782ca8d0dc6c0d94ae619889cf279",
"2f93f0e98517691fb8a30db04fb35c651",
"1207a64a1515928afe34cf4c46e1a2ccd",
"19a858986e34cd9246e55e16cc1faac7",
"2feb403300c71366c14ea17aa63f2686d",
};

std::string _p131_y_str[] = {
"3de01d2bbab8239bcb8fc4a2052c9983",
"b79bb30318d5b37257786b10064c578f",
"28aad11220cb17e0d25dee852bd4030e0",
"77ec9e3a44eea4a5ff3dcd294145d6af",
"14958350b78dec04210a100521bed43e0",
"1545a41126c3ec646de6781d7b12ec833",
"2ef82a07206cef59359c794bd17efef83",
"4159bb2463af3f17c047e2ce3a6aaa96a",
"48862f066652dc01537acd35c2fa4f390",
"b866e9137ccfe38aa6b5d735766e1e8f",
"c88340b9e7e6be462cd4b39f9f8a9192",
"18a442ea1f4066c4e9275f176edd1a07",
"1b636cfc33b6c1e58eab895117085bc9a",
"436f673b28411ad9ebf67c8c68af09cdd",
"40bc7a67bb9fc3d9fd1677445a8e7daf5",
"1c8078a5617bcdd0195cd9650017d5950",
"12e98f7f9e4d85ceff8bf8cabf7225fbb",
"2d1eada5653d18c17e188a8c85d0f96d7",
"3c14ed85739ba24247573e1108ecdaab3",
"2c42cc55042b3176912cd236e03172d5d",
"1eaf2567d575c68e290224001994159b0",
"e3685eb2baae0c51bc43625c185a4b32",
"1efc2edd2b064acebd3d554d1267a1a0d",
"8068821493c5b83f9613c789701ec08e",
"df14396ccf47395abca5352464b42a30",
"af381191209b4fd0f7537416363b7f98",
"3f04a4df1ff385d1c7b0604f018ac55ac",
"37b6ca35cde908e3f58ab9bea4a04ad37",
"1ff0f3fce20eb8dd086aa42b59dbfaf95",
"44bca5522945cf9315221f39a89b43073",
"407746effb5fbd26f4a54509c8769eac",
"2d24257fa2f8e217a03a4befcc2eda825",
};

}

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


std::vector<uint8_t> encode_dp(const DistinguishedPoint& dp, int dbits)
{
  std::vector<uint8_t> buf(ENCODED_DP_SIZE);
  uint8_t* ptr = buf.data();

  // Remove distinguished bits by shifting right then
  // copy into array
  uint131_t x2 = mont::rshift(dp.p.x, dbits);
  
  //int len = ((131 - dbits) + 7) / 8;
  memcpy(ptr + X_OFFSET, &x2, X_TRUNC_LEN);

  // sign bit: 1 byte
  uint8_t sign = is_odd(dp.p.y) ? 1 : 0;
  ptr[SIGN_OFFSET] = sign;

  // exponent
  int len = (131 + 7) / 8;
  memcpy(ptr + A_OFFSET, &dp.a, len);

  // Walk length: 40 bits (5 bytes)
  memcpy(ptr + LEN_OFFSET, &dp.length, 5);

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
  memcpy(&p.x, ptr + X_OFFSET, x_bytes);
  p.x = mont::lshift(p.x, dbits);

  // sign
  uint8_t sign = ptr[SIGN_OFFSET];
 
  // Calculate y component
  p.y = ecc::calc_y(p.x, sign);

  uint131_t a;
  memset(&a, 0, sizeof(a));
  memcpy(&a, ptr + A_OFFSET, a_bytes);

  uint64_t length = 0;
  memcpy(&length, ptr + LEN_OFFSET, 5);

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