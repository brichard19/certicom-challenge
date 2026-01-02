#ifndef _BINARY_ENCODER_H
#define _BINARY_ENCODER_H

#include <cstring>
#include <stdint.h>

class BinaryEncoder {

private:
  static const size_t DEFAULT_INITIAL_SIZE = 4096;

  uint8_t *_buf = nullptr;
  size_t _count = 0;
  size_t _size = 0;

  void resize(size_t new_size)
  {
    uint8_t *tmp = new uint8_t[new_size];

    std::memcpy(tmp, _buf, _count);

    delete[] _buf;
    _buf = tmp;
    _size = new_size;
  }

public:
  BinaryEncoder(size_t initial_size = DEFAULT_INITIAL_SIZE)
  {
    _size = initial_size;
    _buf = new uint8_t[_size];
    _count = 0;
  }

  ~BinaryEncoder() { delete[] _buf; }

  size_t get_size() { return _count; }

  void *get_ptr() { return _buf; }

  template <typename T> void encode(T value) { encode((const void *)&value, sizeof(T)); }

  void encode(const void *ptr, size_t data_size)
  {
    if(_size - _count < data_size) {
      size_t new_size = _size;
      while(new_size - _count < data_size) {
        new_size *= 2;
      }
      resize(new_size);
    }

    std::memcpy(_buf + _count, ptr, data_size);
    _count += data_size;
  }
};

class BinaryDecoder {

private:
  const uint8_t *_buf;
  size_t _size;
  size_t _idx;

public:
  BinaryDecoder(const void *buf, size_t size)
  {
    _buf = (uint8_t *)buf;
    _size = size;
    _idx = 0;
  }

  template <typename T> T decode()
  {
    T value;
    if(_idx + sizeof(T) > _size) {
      throw std::exception();
    }

    memcpy(&value, _buf + _idx, sizeof(T));
    _idx += sizeof(T);

    return value;
  }

  void decode(void *buf, size_t size)
  {
    if(_idx + size > _size) {
      throw std::exception();
    }

    std::memcpy(buf, _buf + _idx, size);
    _idx += size;
  }
};

#endif