#ifndef _P131_CUH
#define _P131_CUH

#include "uint131.cuh"

// ECCp-131 P
__constant__ uint131_t _p131_p = {{0x194c43186b3abc0b, 0x8e1d43f293469e33, 0x4}};

__constant__ uint131_t _p131_k = {{0xe0587d72985b105d, 0xf1fd54b0309e1ab9, 0x7}};

// P - 2, used in the inverse function
__constant__ uint64_t _p131_p_minus_2[] = {0x194c43186b3abc09, 0x8e1d43f293469e33, 0x4};

// r^2 mod P
__constant__ uint131_t _p131_r2 = {{0xa0a3d1806d7fee3d, 0x93a974dd09a9f553, 0x1}};

// 1 * r mod P
__constant__ uint131_t _p131_one = {{0xe6b3bce794c543f5, 0x71e2bc0d6cb961cc, 0x3}};

// ECCp-131 a in montgomery form
__constant__ uint131_t _p131_a = {{0x9893f47bd114555d, 0x2b935e3727d799fa, 0x1}};

// ECCp-131 a in montgomery form
__constant__ uint131_t _p131_b = {{0xe59d7a40bc8648b, 0x55abb9bc04daa824, 0x4}};

#endif