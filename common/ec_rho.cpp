#include <cassert>
#include <stdexcept>

#include "ec_rho.h"
#include "montgomery.h"
#include "binary_encoder.h"

namespace {
// DO NOT MODIFY
// These are deterministically generated and must
// match on all systems.
std::string _p79_a_str[] = {
"1614aeb8e92db814df57",
"4027a11e368e956d92af",
"5167d26cf01aca8bf67c",
"732e3e470ae7638ed05",
"5a1b3cb7c1c5d0631ab6",
"23bb3fd83dbc9352d5e1",
"3b499686e5018ffe3e2b",
"2de349f117f79bf77dbc",
"494bed2ef651d0fe4b6e",
"3a3caadb82edb51b377b",
"391d345ae163b3edb2a5",
"374a4f58940039749c8f",
"60ddf2e5033edc753ef9",
"3663cf9cc6d126c292fb",
"12c308adf7f92a59fe65",
"79084e620ce46337b79",
"5334c520d01bcca5f438",
"271dfc84af3c4191bd53",
"3bbb41fcba037324ef30",
"30a6d7ea510a40bdbd36",
"17fc76572f0e53b9b860",
"627d40ad93ece4b4187",
"291d2f8d459f90aaa910",
"24a5d671bb553afde95a",
"606447988cc6778d59d0",
"11b85dc7f5d8d8fe0db9",
"2e3ddc6562ff7f880357",
"29072a23fad221d674b1",
"3e6099679a21ab541c45",
"430932e0c44731c97bc3",
"30de088da3da65a8da14",
"33c89818c30aa03c38bd",
};

std::string _p79_b_str[] = {
"26ad00a3569c3785004e",
"2320f5a64d245a6e48e1",
"25318ac49ecc98ccb8e8",
"4c9af6b038efb323b91e",
"12bff8347f3f2844d9a9",
"381b2920dfe5941c94c0",
"415ab7e9fe8ba6914814",
"94ee5c57b318e7e416c",
"13dd9410cec72986530b",
"5ef10216635657136ff9",
"ddbf26d833d9c4879bd",
"c9c6249ea691d445b12",
"15557da3f36a91ebfc0d",
"27f8e3c12af04271c534",
"49d1535780bceb4fd6bf",
"1e784e1c9252ac1e790c",
"21215b6fef72e6328835",
"32d94875ecb56b920dcf",
"4ac05cb562e67f279fec",
"43c291909a8afd751912",
"17b829d768e142338e05",
"677d94839904aad3720",
"4fabf6b53722e7c46de3",
"2d4404483feffaddc45c",
"f3031b646c1245e52cc",
"97d2be05cadebd9cb1d",
"5e7195eefb021fbba33b",
"1197fb72ac16016f55e0",
"c14c2e53d858a989c0a",
"169c545411f32d9b958",
"2a0a3c3ec79b5120a83f",
"1110024670d75ffcf4e6",
};

std::string _p79_x_str[] = {
"54899ef7ce3b5d361cb8",
"4bc33690f31706e4e8e5",
"8a7aa1e491dd25a94a5",
"14010418bd9e84a569d9",
"2c4fb3b60398d2d12c53",
"8637cdeb68e39c480fe",
"23f4976bc577229d1808",
"52cbedf28a1080b5a49c",
"8dd23df5b2ae25708fe",
"40de1b1e578e70e105a1",
"34d248f64b8749bb8390",
"539e595a3445fb9ba276",
"2d87f9bf88087d6a4e32",
"1db650dc7b73ff5d5378",
"1182df8955ea22e2def5",
"2d31d8ae49d0e4e6f39d",
"606e79ef78eeae32740e",
"1b418c346b37f7a97a59",
"4ad97e26df47b6433def",
"2d658a8907169c6ad797",
"190a83c768ce5a83a87e",
"16a28c7c65410c877a56",
"20762eae9278d3416bb3",
"3a16b4db55907fcbf4cc",
"425aa926ec569299a970",
"11d66b69557e1315e4ce",
"a11982e2d036e3bea4e",
"29cb7a1bac206bf1deb5",
"ad2b1d3108a65848ae7",
"5edb49b7037cc1fb4c27",
"4a747ec14a28cb917d6e",
"2f1bc0c73a387ec9f353",
};

std::string _p79_y_str[] = {
"3c1b2c03272abb046764",
"5e09e34a1d096154e348",
"3622ba45d05f3883ad96",
"46726cf1cbb64ded7b2",
"616b68190100da4deda6",
"2eebdc9c8f1d396905e8",
"112af321f33bdcbd2fb4",
"3b70a6ef3c1040e40af0",
"4d93ff2237e9aa917799",
"41d6422689439cd1b94c",
"52bb45b63667cfdfd28e",
"23ef519409069fe9b0d",
"155a1fd14713bf1d2f5f",
"872fb24cd95f034faa5",
"6ec7ac6d44b16773b87",
"3184a893a3f118ddaa3f",
"37811ffe9991a00e898a",
"2ffe839e5e064179f297",
"1a301f19cdd505df601",
"16d1686056f7bd1cbf3a",
"6266ec41cd8d7a0ebece",
"48f6bb6d28762cd255ac",
"9e7812f58e66d258fcc",
"1143c72c689382ed0dce",
"490151fdf28107e90860",
"50b71ceacd964bc208bd",
"368b835d236b6ed394e4",
"5c3479c1793fd1bb7573",
"14f915988833b96c0d73",
"2dcd409aa825ccb09e45",
"2bc38b1cf02e0c431a71",
"6244f6aedfa3b2607674",
};


// DO NOT MODIFY
// These are deterministically generated and must
// match on all systems.
std::string _p131_a_str[] = {
"142b5e364eb43a80dd5f2e30a9e501b7d",
"1fe4bba7faba3f6a5fba9a1630e23bac5",
"16aa058d69f60503d5cd158ba16a76b0d",
"206113ff3eb9914d541b96ab112ea4c55",
"1b7802fb7705aa46d48a12d9e7a7f7e9d",
"14bf17484621d5f050e3ac880e5fe21e5",
"22a3cd2359345a49dc22299e52381562d",
"10372ac9e6da3d535ce2f82f15b283b75",
"1c060551e9bf450cd92fe80074395f1bd",
"18c2c438e635f7765b17b50af66719905",
"ca207789acf9a8fd9ab990bf74e6514d",
"2a2914ddfbf434595f6a68cedac233a95",
"1a5ecbcada7a8ad2dc3854d4349db74dd",
"1a34c276954023fc5df74c82010c62025",
"1cbacd2e35c145d5da144f0c0dd1e5c6d",
"14ee370057b0f65f58c8b5c7b492349b5",
"19a813ad3b90fb98d7f4c8a20619807fd",
"1949fc9f6449db8259ddc6f386a43b745",
"150aee1d1fc2dc1bd1720cf09b271778d",
"167186cf5a7a03655b293079c797068d5",
"a2652261e1c175ed3fce712de313ab1d",
"1a66b60f1b1c9e08501d524c3fc325e65",
"135445de9d4ddd61d5e4ffed4cf27a2ad",
"1eff4f374b78db6b589d8fb32885297f5",
"187f222911f55e24dffe997aeaa965e3d",
"1cd58bab9841eb8e50082c32653da1585",
"22a9e904a19bc9a7d9361fee99188ddcd",
"178dd396b896fe7150c08c0afc511d715",
"1532570085b64fead7190d38b1868215d",
"6779f6131034414596f36e8d0282dca5",
"1716c3022ea620edd14d9f80ddbdd28ed",
"12d86992d67dec7759537f4c982f62635",
};
std::string _p131_b_str[] = {
"1ae12a08b93639bf96438f14dfe833301",
"177387c9783cda411552e421879d836c9",
"571d6004137e67298ade8c37577c5b91",
"1b2418bb680dc4541087799f31d30c159",
"e23cda189bf99e5911aadb8d2d7a8821",
"1f9a1230cf814d27102e560ba1422cfe9",
"1204edfb4cd1841890342921fe2b6b8b1",
"13276af2d891a4ba1dbdab32a9d076279",
"1c6e2902c11dd50b91ffcc8c8b5a9ed41",
"2565fee1bb64fb0d12d37dde19a777909",
"14dbb0326d00bcbe9d632c298510d25d1",
"1fb878b4f14d8020106b44fac234c1399",
"2132d577b9826b3192589eb895bd96261",
"185f3fa227c963f3106ef8adc129e3229",
"11a1289f455710649199caab70947a2f1",
"1be7251af382d686124a5cea097c6d4b9",
"2099836ff7dedc579151ae117a8d0e781",
"14b23ad1435007d9188f17162c65efb49",
"23bfe2a4d425ff0a927bf2dbb362e3011",
"16a8e532933327ec1888e2526263fa5d9",
"209415f58be4a87d9234a508de9587ca1",
"2481ea89df5a66bf151af9eed4381d469",
"181d7d5dd27f08b091a1a449ec688cd31",
"1343eaf8571ff4521ca295a223e7e86f9",
"e7f35d370054fa3968e6eaea2e3821c1",
"e695974d00a00a51e26171035bcebd89",
"1cf4439b1433ad56968b25ea86d1f7a51",
"2a29ba01f8cabbb81d0d4c313944b7819",
"14fe746be97251c996a0f4c204c37d6e1",
"1897954b4c40558b1c3cee28f250db6a9",
"1950df91e7d56cfc9ab83e2fda0ba3771",
"18423dd8c474fe1e11f6481e41f6e7939",
};
std::string _p131_x_str[] = {
"4d29b7324d05db8ecef3cbc431c57b9a",
"135d3712059496347e42ca1bdc1430135",
"276e501d8b98bb114e20498079063acda",
"33040bc01a4e93f12789f4c40737577f",
"3d85797172ff14274b3db1a8b96da5828",
"42f1b1259c5c66a8d23f198ff802fe561",
"42bc0f3e59a5d9e59f9479b928dc76c99",
"162162788aa781c8157b1c8e169be5a17",
"427f0a96d50c7bc9de795067c27b0e82a",
"1c478dc55e0632bb70d369533b8c31fbb",
"3a927d4ca9d6b5e67543058be3f15e78f",
"62ecbb76630462c8850a798dcf1e5c46",
"24d8f23cf1e6b8bf084d6008d93a3f5",
"2f4681411a9ed48b30f273c3763087f40",
"bf8dd6d1f7e30347188f2373f5048e90",
"2abfc7085b06ad53d04c45f5b938f467d",
"2d07bee5b7621660e34e51f780dc293ae",
"4091a7ccc939bb1c8353fe00213b4f172",
"3cc744f6e889f051862a7c1425430b1f",
"2f8a6f88c7ace9fa4acc8fede6ec3dc2c",
"27850bafa1b838bb943e95af0df30d7e1",
"4599da7d0064f98d358c61eec8240f44a",
"43f52802b539f00b9dfcdb9f56b4ae528",
"32458d364e5af94f7bd0c0d6aefecdb32",
"1cec55d33bd3ab5dcac973118ed37fa9",
"472fd48e9fae3b1e8d96091bce8efa5c1",
"e65e4d17b1c05a379d185ad5f4c8a7a6",
"4dfa0ad8553ce4728dcc509e40e8a4b7",
"44a1ccea2cd4f76d8fb42f64c0a4f9277",
"5b1be0dd20ef16024b53cb6d38a63be4",
"2f8300531ca303305baf2c28f11c8b572",
"1337af6c71437663fe3244d6da3e7cae3",
};
std::string _p131_y_str[] = {
"105af2008563de90bba00cc6048953294",
"10f5c7c1a77dd9a00a1ddbcdf962fc57f",
"f8317bf20e3161be423e036e967933d8",
"21e6f9b05be250db3b13ff67a5f152102",
"45cd26bef9b2986693f3954109c3c523c",
"142731555df1e2085fa57fc406f23a9a6",
"5002e06750b6974b352d8f9a1b821c75",
"2e0938e022c022f6f1cded6fed6577876",
"34eb01e413e0c4ddd9252cc7b3022e712",
"308536d7315e77aec8a53745544dd97a2",
"1beaa6db6ee895248bb5eddb08497482b",
"2ca08734dabe3ce0771fb9bccc9807d4a",
"79650f6a2ac37beea155b4451c1e00b4",
"13345f72865473c457494cfe9c84f2bd9",
"4578b35807f820d96f3bb8a55c742795f",
"17e13da6e5cb877bc30e8f3a45dffdf7c",
"37219a6f1a473491e1bc0b89d0d513c0",
"129f728011768536e910adc8c8df91d67",
"8f8bb52ff5d667a78c57e77af671c6df",
"b05556fc59fcf134571207b6980c473b",
"202a46c5ab5c5c9107fb4ee68d77fc166",
"10c24995b7982090bac4cf1fd6e4f6d94",
"34208cd9044bb7dd590ffb08e3a6dca16",
"43a6cd235694b9c2adb4203ccf5fc118d",
"169c27fd7eea9ad114a69627d39e9a844",
"111934e822b0c2ef647b591815dcd879c",
"3ea80981d2cd4e0d9c3d9bcdfda37c656",
"12c65b5b295928d6b93078ccf66b9c77f",
"41180f2dd884bf79b86f037f6a08b748f",
"bc2ae81c926b5197f1f2dfde4ae9c243",
"44daeb45c4bb670cf3dec9048df711659",
"3ccb58d478406ebb23c85770ea889caed",
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