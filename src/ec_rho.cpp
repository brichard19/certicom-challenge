#include <cassert>
#include <stdexcept>

#include "ec_rho.h"
#include "montgomery.h"
#include "binary_encoder.h"

namespace {
// DO NOT MODIFY
std::string _p79_a_str[] = {
"43f08e2cbaee90bf2d56",
"423768de2b386bd6ff5a",
"5a7ee45173b1266a956f",
"2b6b3c4004c46dbf56f0",
"614fce830be06ded1db0",
"5d7f63c4ff45d9442007",
"2d9953723be1b6a60118",
"2319a336e77f25f7c245",
"301f9b52ce930e6f102d",
"112e9c800b22c52c681d",
"142d901f2bb9474e4e66",
"402aff7bbf8f1dd016b8",
"deb9fb9e665a35de168",
"45db7733b44863ae710d",
"452471abe3b506828108",
"5d32a439ddb48bd8b1ec",
"2a35a88c3dc3fd1f9a3d",
"276f6fefafd2acbb1f92",
"ddef24ebf51a451b958",
"465762450665256da72e",
"c48f0ebbfb20a79071e",
"51e6c7a7852818b4d622",
"2a91ed1305dd950a7aa7",
"b1bd007b4ba70d225cb",
"4f3cb35c077e9a272383",
"3456fe90294c5a3833ae",
"5751c265465321340d32",
"278a31ff7844edf1eb51",
"1f2cbbe912ca5812f402",
"33fe34c2c6324b0e7e42",
"57dc28d8df2c92ef2dd5",
"10701fb26dad6b6520a9",
};

std::string _p79_b_str[] = {
"4e229b6aa4c02c8b5cbb",
"1381c734c0fda30cd865",
"4a1f398aee869afb62fa",
"52b649d4cbfcc232a3f0",
"421c8640ec1d56847f77",
"357c70b972a9bbd7460f",
"5b5d20e6db69ebf32a1f",
"4061089da95f6a2f0e13",
"3c552875606633b0b591",
"2b4774e63ed4da9796a4",
"61f84f2ee3ef5a177996",
"1713c414ff8a3f5b9cb7",
"3be61094b2aac492c2ae",
"567b4cdc872e2058188c",
"1aabc3a72d50671e6f22",
"21c97de355ed129aa377",
"431febe00f76582ff5a6",
"4d02aff3c06e5b7f0517",
"1271811bcac2b48a5885",
"c60bbdd178783ed43cc",
"2849eeaa759f607b7289",
"59c1e9e8c19199fdf94b",
"8b47a554b518e9c3cfc",
"2e54ffe4bfdef5847ab3",
"52764f05998fe9338977",
"37548825b5a6ceed7194",
"27443884674c3e4f31f0",
"5880dd804d35777aee26",
"40839437be868e40d221",
"53686d719a5ccd535278",
"2230a7da912828c90f97",
"32f7829c00c89535a963",
};

std::string _p79_x_str[] = {
"a69ce680f1dfbca52ea",
"7a0a4cee13b7b90105d",
"2b25d825898a426e4e4d",
"24ec1a3242da90df4284",
"36a343b07dd8fd54ea8f",
"10a807a85a67b534e84a",
"40e28bded898b8587277",
"1ab7fcb2b9713ee2e88e",
"4027f8f55e859470f2ac",
"452c4bf742a0ce71ee3c",
"18a64a7778d9ae3e34a6",
"5ef339c8478f4ff1f246",
"593b46e974a44f019553",
"59e56d029e23f622ea71",
"5606c8bc479e93b666aa",
"114179d753fc0cf387a8",
"6645e03cc1e79625aea",
"5c13bc2f0827caeb3b89",
"d5a52dc9dfc987dc8e9",
"305a503259229f5e77fa",
"5eb8e3c485a70a167644",
"234f97392bf451836a96",
"2dafb3ca7676472e995b",
"181c62e3bb7374faa4e6",
"c955d2ba08799c88c21",
"1d2a6cc95a4fa2d76899",
"56c8c5513d7c8dc6ab2",
"408f4a578b4447164772",
"36d87daf5d3d012d381e",
"45bc5d971de67765c07d",
"688c8ff68283ddfaa4f",
"2f3d872f77eff543d590",
};

std::string _p79_y_str[] = {
"190e3f4f89c7145c5e91",
"b0b493e9f63d8d9068c",
"2d414b2ca81fab6d26c9",
"32da810fc2fec798c99a",
"127225b332f8fe53c65",
"16bb2baffe89c6b1d732",
"fdae995f0bcbfc4d5ae",
"240c682b800eeef693db",
"fbb96bd6a384a4324b9",
"2cd9d647611ffba9183d",
"308446639a8d7fd86c1f",
"3de4e46436cdf27a0858",
"5a4338c5c361c012fe0e",
"161af31abbc1e7b9c6d1",
"3c7e5b4adf76defe1ea5",
"220731bf3f22a496565b",
"134d7ee43a079a1ed6cf",
"24f0bea3dae265a392f7",
"41fbb1d1a83436648185",
"2fa11ab3610eeccc7f4c",
"92c6ba579e8e972e299",
"58504c45517a54231bf3",
"2bad2795f359eecb37af",
"35396c3a1f8549337611",
"62644e29d3a9827e73cf",
"37e325b6c63fe3d2a405",
"1b45e552aa37ee44cbc5",
"5e86e98c3bffdb39909e",
"b4b3ae8c00a47aacee7",
"55e15abc37a3f47ad6d5",
"4e1d2ec680fcc32dc79d",
"51a2d33ccc8a574b9440",
};

// DO NOT MODIFY
std::string _p131_a_str[] = {
"45b422da3150d6baa3212a4fbd777c8dc",
"43d7d063a1ad094462c115c4367a72c1c",
"2b6a72e3311565d7b5ca16924b8ac0daa",
"3d28664312f876b9b8a4c59d01c9190c6",
"2a634e204490887d07f4a312c01a8f76e",
"34cede5db4a5d8b4b81c04ad4b255e10e",
"44313cb192d87b182bea8707df84fdf0e",
"12d2156f37db1da4cc39ce9d97503385c",
"3b4348fe24bd3fd61ec76f69214573cb0",
"9e7424bf97c92ad978a61a4c3e0729d2",
"7431797982d1ffd4e04c3cb95a698ae2",
"29af73e3d9812a71b719fcf33bff0fb8a",
"23b318463d680fc11804117162be9020e",
"37ee473d9801d4c81fc995112c0655df6",
"2865fb3b65de402fed45b8a03543db42a",
"485c69473bd4b6488ff9c4fd116168126",
"1fc6ebfb1273f8de0af074ca0d37e2f78",
"2f7928813311144bef626f4f6d5301e47",
"271dc9193ee76f8befa46e5c85f72332e",
"efc25830af8cb81a280de8a73a79acec",
"14b7d531c5ddf5503ad9c2a6a39345eac",
"3ad6d1783e7c0bb7d56b9e28c50ffb124",
"47adaff61b9ec675ccf619d9233fc5db2",
"3fe1e8ae7a49c72210b330c2c863c0dbc",
"2c67275049f0e7d78c6181b4278a42069",
"1a25c1afe20af178a9e8b9c1e183704a4",
"3dd6f4d1c20da81547818d1d6aa0a6d42",
"46797fd3f7b26e35d3fa15098ae40af98",
"230963329a0396de63622b51575204bec",
"337715abc171332ec41d9055b26f5ada8",
"39cec2ef6103fb31532295ed2379dd43c",
"39f7ed0bd1961cd5d46146abafdcc621b",
};

std::string _p131_b_str[] = {
"2474250a6f792807345488df841cb61c3",
"3820da3eedf2582025eada5128f6adf90",
"20f56f2e751258d03271a39c9bb26527",
"401bade6c981bbb6a59901f0342dfe550",
"3aa227658aef8a21b7fe0ae27eb4a39a1",
"3746e2ee8afdf4fe510e125e2433b4f3",
"37ffb5d8c1e7752ea7f006f22c06ba82",
"2420e328d8cd15b5e342c277a44e6292a",
"150b5d7568090d8f49681413d5b57a80d",
"242843db808e47a2f232963a68dafc410",
"2f07f519d69620f42072341d14550f1a4",
"4858cf84918717368c11dcc3ab2dbd5f0",
"3799037c5039a8d513298ee3dedcd1a94",
"3395ad053abf9425be9dc0a3d19e9f07",
"6d092debd95aefe58c784dd2386ecfc0",
"1dab7d35b0e0b905c1bea6b1e8392e83e",
"1b9d62e9fc51e42f2d72d7f2a4f5127fe",
"e50e104a488141e9e175e020e8a095f4",
"2cc1d33d7acd8f0a34f52ab2b57864809",
"77a3593266d510fd4df64dec7292b57f",
"2063f1a45ebbdf315e3a266f2b413723e",
"2b9a427c0e7c5c28e0712c6e46f2ba269",
"44ad72d1689cccb72e084b87d96052682",
"1c5901ad2ee135bdb67f1826efc5b1954",
"2803625d5bc8a0791a1ebea34050d5988",
"16625db2205a84f99ae6a550729ece6ba",
"16bc9acddd6a385bdf7ed2be4d6009a93",
"1749e0d9a44c65e8f0bc0eeadd99eabd9",
"8396993337d86842945edbfc6abe0aeb",
"396625da14d4882be75480b82428d42bd",
"32d8d34a8797b0ea5731f19a3f2cea0ff",
"1e6e926eeb847929d8d57df41bfb87856",
};

std::string _p131_x_str[] = {
"26504caad9a173e80e53954fcb8f7d763",
"48dfe3a87c30439b8884f0b104e0e4598",
"32b19271d54efd3142f167c6edf1eeced",
"c3904474213a3af8a9462f64923856af",
"19740500486a55ca2b57595faa0f360de",
"3e44e55e6af126d3f9ff31f5d03526440",
"10819a0bc8cd81a63b9814931a9adaae9",
"700ee65a8a9ddd407bc6f1c7a9e6cc72",
"283c024cd43d12311435cbd3997246917",
"7293b406f32ddb1bd22ca3d667e7cf5b",
"3ee962c9011846d1962e0deeb3c052962",
"459db76f1309136125a182c2474611f58",
"87e1d51aaf6c0a097604fcf776459cfe",
"3a48173a942f3752db0baa1d5b9ac7e16",
"43446aae5b008d38f6eff4d5e5bfebb76",
"ed47a07c80272a061c681a89bddabc70",
"381d89ac839a31348e98e391c0b3a89d",
"1ece27b3005fcb98a1453cc534c5682f2",
"3089060a3d19ad9e5976c77fb5753d8a",
"43c050f4c11d3e9bd1d2be481c32a0bbd",
"418418122c1d6d66063a25f846d8a6731",
"470930803ceedeb983824fe2dc284ff16",
"75697b152af06342163998e5ad42c9a9",
"1adccc5ed5ecfe0d50110db8ad9960aa5",
"aacfeade9afda642c3a1c2286de0b0dc",
"2a8b0e83adc660886b27f4612b59a8521",
"235941880b53e8b3502c8092674f1b83b",
"470e5dad7dc10c4bb2722cc118e99def3",
"2791399d90e3bdc7206fddcbdb7d65a52",
"15f7ec63fd5c55cc743aad4453c0abece",
"fb4d39a4c1e1d3e0a5ca31607f0f8d80",
"3c7a99b4251551d467a176dd8b760c971",
};

std::string _p131_y_str[] = {
"3555082d885f9c405b6be6b4fe2df905c",
"2fa5c76c543be6533c59d78d49440cbd0",
"44db48e5cb3abbc2ffc30bc5fdae8ba8",
"25f780836c3db01ed5e6a40183d82c5ec",
"1f9ded38e77a01d456ec442c8fadf1b75",
"3104002397d923d48c133140e46bca8c9",
"3f5670a9316031e2e011dace5fa60ff6d",
"3d0cecc9d8309998bb531c33b187970c",
"30386666fb0149d349239e4fda746b466",
"18b2e71c4c10279637b2a9e62763d260c",
"2a8cb5d6d1ab4fd2efc097510ea06e8e0",
"4148ae057d95fd7b69f1d609cc2183fe8",
"4611a2ac562da59406ce1d586e31786e3",
"142f240e189b5365148d53a2f2313ae67",
"ee87c3081179cb03eebd4d088ca7d0da",
"23f64bbc0a1f3d74cb1928f569de1666f",
"3621ae01dd4448601da2a93175e110699",
"cb40ecd8ce9f43a5ea3e9224446c48f2",
"430effc41c28cf45234349aa149b00bf2",
"28e6ad38ffb0216803abd657417f3e237",
"8510a56c326e55a9f1a806d00d316645",
"26017715f60345bdc72c289000c0bad0b",
"475d03d0136ff80d99352a50ddfe2b4be",
"274d8e87fb2ef4509a399acdd4a5c9814",
"361aad4c33e869f0969c877730796e967",
"25aeef8e4dbedb5a8f015196932f37a19",
"3f8aaf5182ef7c4c6ef9dc2706ae0e626",
"d29707b5072fb26262d119032d9a07ce",
"e4fe2db5bce8331dc89b61140137a1d2",
"15a19c9abc070efce1ec3eb8b3ac30d3",
"3f525ffa1dbd3604ee4117ca41a966cbc",
"31f0fd9355b0462b6f3ffc6d71ea4b499",
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