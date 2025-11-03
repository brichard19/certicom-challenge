#ifndef _GIGALIST_H
#define _GIGALIST_H

#include <fstream>
#include <string>
#include <stdint.h>
#include <filesystem>
#include <stdexcept>
#include <vector>

template<typename T> class GigaList {

private:

    std::string _full_path;
    uint64_t _count = 0;

    std::fstream _f;

    uint64_t append_internal(const T& value)
    {
        _f.seekp(0, std::ios::end);
        _f.write((char*)&value, sizeof(value));

        uint64_t idx = _count;
        _count++;
        
        return idx;
    }

    std::vector<uint64_t> append_internal(const std::vector<T>& values)
    {
        _f.seekp(0, std::ios::end);

        _f.write((char*)values.data(), sizeof(T) * values.size());

        std::vector<uint64_t> idx(values.size());
        for(int i = 0; i < values.size(); i++) {
            idx[i] = _count + i;
        }

        _count += values.size();

        return idx;
    }

    T get_internal(uint64_t idx)
    {
        if(idx >= _count) {
            throw std::out_of_range("Index out of range");
        }

        _f.seekg(idx * sizeof(T));
        T value;
        _f.read((char*)&value, sizeof(T));

        return value;
    }

    void init()
    {
        std::filesystem::path file_path = _full_path;

        auto parent_path = file_path.parent_path();

        if(parent_path.empty() == false) {
            std::filesystem::create_directories(parent_path);
        }

        if(std::filesystem::exists(_full_path)) {

            // Verify the size is a multiple of the data type
            uint64_t size = std::filesystem::file_size(_full_path);

            if(size > 0 && size % sizeof(T) != 0) {
                throw std::runtime_error("List size not a multiple of data type size");
            }

            // Open for reading and writing
            _f = std::fstream(_full_path, std::ios::binary | std::ios::in | std::ios::out);

            _count = size / sizeof(T);
        } else {
            
            // Create and open for reading and writing
            _f = std::fstream(_full_path, std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
        }

    }

    void sync_to_disk()
    {
        // Close and re-open
        _f.close();
        _f = std::fstream(_full_path, std::ios::binary | std::ios::in | std::ios::out);
    }

public:
    GigaList(const std::string& file_path)
    {
        _full_path = file_path;

        init();
    }

    ~GigaList()
    {
        _f.close();
    }

    uint64_t size()
    {
        return _count;
    }

    uint64_t append(const T& value)
    {
        return append_internal(value);
    }

    std::vector<uint64_t> append(const std::vector<T>& values)
    {
        return append_internal(values);
    }

    T get(uint64_t idx)
    {
        return get_internal(idx);
    }

    void sync()
    {
        sync_to_disk();
    }
};

#endif