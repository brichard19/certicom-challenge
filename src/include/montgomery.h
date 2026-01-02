#ifndef _MONTGOMERY_MATH_H
#define _MONTGOMERY_MATH_H

#include "uint131.h"

namespace mont {

uint131_t to(uint131_t a);
uint131_t from(uint131_t a);
uint131_t add(uint131_t a, uint131_t b);
uint131_t neg(uint131_t a);
uint131_t sub(uint131_t a, uint131_t b);
uint131_t mul(uint131_t a, uint131_t b);
uint131_t square(uint131_t a);
uint131_t inv(uint131_t a);
uint131_t sqrt(uint131_t x);
uint131_t rshift(uint131_t x, int n);
uint131_t lshift(uint131_t x, int n);
bool less_than(uint131_t x, uint131_t y);

}; // namespace mont

#endif