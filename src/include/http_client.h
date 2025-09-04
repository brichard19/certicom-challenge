#ifndef _HTTP_CLIENT_H
#define _HTTP_CLIENT_H

#include <stdexcept>

#include "json11.hpp"

class NetworkException : public std::runtime_error {

private:
  std::string _msg;

public:

  NetworkException(const std::string& msg) : std::runtime_error(msg) {}
  NetworkException(const char* msg) : std::runtime_error(msg) {}

  std::string msg()
  {
    return _msg;
  }
};

class HTTPClient {

  int _port = 80;
  std::string _hostname = "127.0.0.1";
  std::string _url;

  std::string post_json(const std::string path, json11::Json& json);

public:
  HTTPClient(const std::string hostname, int port);

  int submit_points(const std::string&curve, const std::string& encoded_results);
};


#endif
