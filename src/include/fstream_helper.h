#pragma once

#include <cerrno>
#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace ifstream_call_detail {

[[noreturn]] inline void throw_error(std::istream& stream,
                                     const char* command,
                                     const char* source_file,
                                     int source_line,
                                     int error_number)
{
  const auto state = stream.rdstate();
  std::ostringstream ss;

  ss << "File stream command failed\n"
     << "  source: " << source_file << ':' << source_line << '\n'
     << "  command: " << command << '\n'
     << "  stream state:";

  if(state == std::ios::goodbit) {
    ss << " goodbit";
  } else {
    if(state & std::ios::badbit)
      ss << " badbit (unrecoverable I/O error)";
    if(state & std::ios::failbit)
      ss << " failbit (operation failed)";
    if(state & std::ios::eofbit)
      ss << " eofbit (end of file reached)";
  }

  ss << " [rdstate=0x" << std::hex << static_cast<unsigned int>(state) << std::dec << "]\n"
     << "  characters read by last operation: " << stream.gcount();

  if(error_number != 0) {
    ss << '\n'
       << "  system error: " << error_number << " ("
       << std::error_code(error_number, std::generic_category()).message() << ')';
  }

  throw std::runtime_error(ss.str());
}

} // namespace ifstream_call_detail

namespace ofstream_call_detail {

[[noreturn]] inline void throw_error(std::ostream& stream,
                                     const char* command,
                                     const char* source_file,
                                     int source_line,
                                     int error_number)
{
  const auto state = stream.rdstate();
  std::ostringstream ss;

  ss << "Output file stream command failed\n"
     << "  source: " << source_file << ':' << source_line << '\n'
     << "  command: " << command << '\n'
     << "  stream state:";

  if(state == std::ios::goodbit) {
    ss << " goodbit";
  } else {
    if(state & std::ios::badbit)
      ss << " badbit (unrecoverable I/O error)";
    if(state & std::ios::failbit)
      ss << " failbit (output operation failed)";
    if(state & std::ios::eofbit)
      ss << " eofbit (end-of-file indicator set)";
  }

  ss << " [rdstate=0x" << std::hex << static_cast<unsigned int>(state) << std::dec << ']';

  if(error_number != 0) {
    ss << '\n'
       << "  system error: " << error_number << " ("
       << std::error_code(error_number, std::generic_category()).message() << ')';
  }

  throw std::runtime_error(ss.str());
}

} // namespace ofstream_call_detail

#define IFSTREAM_CALL(command)                                                                     \
  do {                                                                                             \
    errno = 0;                                                                                     \
    auto& ifstream_call_stream = (command);                                                        \
    const int ifstream_call_errno = errno;                                                         \
    if(!ifstream_call_stream) {                                                                    \
      ifstream_call_detail::throw_error(ifstream_call_stream,                                      \
                                        #command,                                                   \
                                        __FILE__,                                                   \
                                        __LINE__,                                                   \
                                        ifstream_call_errno);                                       \
    }                                                                                              \
  } while(false)

#define OFSTREAM_CALL(command)                                                                     \
  do {                                                                                             \
    errno = 0;                                                                                     \
    auto& ofstream_call_stream = (command);                                                        \
    const int ofstream_call_errno = errno;                                                         \
    if(!ofstream_call_stream) {                                                                    \
      ofstream_call_detail::throw_error(ofstream_call_stream,                                      \
                                        #command,                                                   \
                                        __FILE__,                                                   \
                                        __LINE__,                                                   \
                                        ofstream_call_errno);                                       \
    }                                                                                              \
  } while(false)
