#!/bin/bash

find . \
  -path ./third_party -prune -o \
  \( -name "*.cpp" -o -name "*.cc" -o -name "*.cxx" -o -name "*.h" -o -name "*.hpp" -o -name "*.cu" -o -name "*.cuh" \) \
  -print -exec clang-format -i {} +

