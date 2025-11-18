#ifndef _ECC_INTERNAL_H
#define _ECC_INTERNAL_H

#include <string>

#include "uint131.h"

struct CurveParameters {
    uint131_t p;
    uint131_t a;
    uint131_t b;
    uint131_t n;
    uint131_t gx;
    uint131_t gy;
    uint131_t qx;
    uint131_t qy;

    uint131_t k;
    uint131_t one;
    uint131_t two;
    uint131_t p_minus_2;
    uint131_t sqrt;
    uint131_t r;
    uint131_t r2;
    int bits;
    int words;
    std::string name;
};

extern CurveParameters _params;


#endif