#include <exception>

#include "curl/curl.h"
#include "fmt/format.h"
#include "http_client.h"
#include "util.h"

extern "C" {
#define CURL_STATICLIB
}

namespace {

size_t write_callback(char* contents, size_t size, size_t nmemb, void* userp)
{
  ((std::string*)userp)->append((char*)contents, size * nmemb);
  return size * nmemb;
}

}

std::string HTTPClient::post_json(const std::string path, json11::Json& json)
{
  std::string json_str = json.dump();

  CURL* curl;
  CURLcode res = CURLE_OK;

  curl = curl_easy_init();
  if(!curl) {
    throw std::runtime_error("Error initializing curl");
  }

  struct curl_slist* slist1 = nullptr;
  slist1 = curl_slist_append(slist1, "Content-Type: application/json");

  std::string full_url = _url + path;
  std::string response;
  curl_easy_setopt(curl, CURLOPT_URL, full_url.c_str());
  curl_easy_setopt(curl, CURLOPT_HTTPHEADER, slist1);
  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
  curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
  curl_easy_setopt(curl, CURLOPT_POSTFIELDS, json_str.c_str());
  curl_easy_setopt(curl, CURLOPT_VERBOSE, 0);
  res = curl_easy_perform(curl);
  curl_easy_cleanup(curl);

  if(res != CURLE_OK) {
    std::string err_msg = fmt::format("Connection to {} on port {} failed: {}", _hostname, _port, curl_easy_strerror(res));
    throw NetworkException(err_msg);
  }

  return response;
}

HTTPClient::HTTPClient(const std::string hostname, int port)
  : _hostname(hostname),
  _port(port)
{
  _url = "http://" + _hostname + ":" + std::to_string(_port);
}



int HTTPClient::submit_points(const std::string& encoded_points)
{
  json11::Json req = json11::Json::object{
    {"points", encoded_points},
  };

  std::string response = post_json("/submit_points", req);


  return 0;
}
