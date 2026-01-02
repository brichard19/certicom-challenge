#ifndef _LOG_H
#define _LOG_H

#include "fmt/format.h"
#include <iostream>

template <typename... T> void LOG(fmt::format_string<T...> fmt, T &&...args)
{
  std::cout << fmt::format(fmt, std::forward<T>(args)...) << std::endl;
}

#endif