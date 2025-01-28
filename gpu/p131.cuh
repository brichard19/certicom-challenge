#ifndef _P131_CUH
#define _P131_CUH

#include "uint131.cuh"

// ECCp-131 P
__constant__ uint131_t _p131_p = {{0x194c43186b3abc0b, 0x8e1d43f293469e33, 0x4}};

__constant__ uint131_t _p131_k = {{0xe0587d72985b105d, 0xf1fd54b0309e1ab9, 0x7cfd70cf}};

// P - 2, used in the inverse function
__constant__ uint64_t _p131_p_minus_2[] = {0x194c43186b3abc09, 0x8e1d43f293469e33, 0x4};

// r = inv(2^160) mod p
//__constant__ uint131_t _rinv = {{0xea2b4033, 0xc3523a9c, 0xabc154a6, 0x3958a1b5, 0x02}};

// r^2 mod P
__constant__ uint131_t _p131_r2 = {{0xf95d709f92600513, 0xf3d6fa1fb65ef639, 0x3}};

// 1 * r mod P
__constant__ uint131_t _p131_one = {{0x6e7743da32b6d0c7, 0x88c614d64c1a8f0b, 0x00}};

// ECCp-131 a in montgomery form
__constant__ uint131_t _p131_a = {{0xe7f7f250cee8709a, 0xacd15fe1a8ec1522, 0x00}};

// ECCp-131 a in montgomery form
__constant__ uint131_t _p131_b = {{0xc85087e5ab4eca9e, 0xde124657d7ba5851, 0x2}};


#endif