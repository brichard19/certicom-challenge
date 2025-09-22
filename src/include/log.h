#ifndef _LOG_H
#define _LOG_H

#include <format>
#include <iostream>

template<typename ...T>
void LOG(std::format_string<T...> fmt, T&&...args){
 std::cout << std::format(fmt, std::forward<T>(args)...) << std::endl;
}


#endif