#ifndef _GIGASET_H
#define _GIGASET_H

#include <vector>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <cassert>

template<typename T, auto CompareFunc, auto HashFunc> class GigaSet {

private:

    // Number of buckets is hard-coded
    static const int BUCKETS = 16;
    static const int BUCKET_INIT_SIZE = 256;
    const size_t BASE_OFFSET = sizeof(int) * BUCKETS;

    std::string _full_path;
    FILE* _f = nullptr;

    int _bucket_size = BUCKET_INIT_SIZE;
    std::vector<int> _counts = std::vector<int>(BUCKETS);

    static bool compare_less(const T&a, const T&b)
    {
        return CompareFunc(a, b) < 0;
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
            uint64_t size_bytes = std::filesystem::file_size(_full_path);

            size_bytes -= sizeof(int) * BUCKETS;

            if(size_bytes > 0 && size_bytes % sizeof(T) != 0) {
                throw std::runtime_error("Set size not a multiple of data type size");
            }

            // Open for reading and writing
            _f = fopen(_full_path.c_str(), "rb+");
            // Calculate current bucket size
            _bucket_size = size_bytes / sizeof(T) / BUCKETS;

            // Read current counts
            assert(fseek(_f, 0, SEEK_SET) == 0);
            assert(fread(_counts.data(), sizeof(_counts[0]), BUCKETS, _f) == BUCKETS);
            fclose(_f);
        } else {
            
            // Create and open for reading and writing
            _f = fopen(_full_path.c_str(), "wb+");
            // Reserve space for bucket counts
            assert(fseek(_f, 0, SEEK_SET) == 0);
            assert(fwrite(_counts.data(), sizeof(_counts[0]), BUCKETS, _f) == BUCKETS);
            fclose(_f);

            std::filesystem::resize_file(_full_path, BASE_OFFSET + BUCKETS * _bucket_size * sizeof(T));
        }
    }

    // Merges two sorted vectors into dst
    void merge(std::vector<T>& dst, const std::vector<T>& src1, const std::vector<T>& src2, std::vector<std::pair<T, T>>& duplicates)
    {
        int src1_iter = 0;
        int src2_iter = 0;
        int dst_iter = 0;

        int num_duplicates = 0;

        while(src1_iter < src1.size() && src2_iter < src2.size()) {

            int cmp = CompareFunc(src1[src1_iter], src2[src2_iter]);

            if(cmp < 0) {
                dst[dst_iter] = src1[src1_iter];
                src1_iter++;
                dst_iter++;
            } else if(cmp > 0) {
                dst[dst_iter] = src2[src2_iter];
                src2_iter++;
                dst_iter++;
            } else {
                duplicates.push_back(std::pair<T,T>(src1[src1_iter], src2[src2_iter]));

                num_duplicates++;
                // Keep the one that is already in the list
                dst[dst_iter] = src1[src1_iter];
                dst_iter++;
                src1_iter++;

                // Skip this element
                src2_iter++;
            }
        }
        
        // Either one or both iterators are at the end
        assert(!(src1_iter != src1.size() && src2_iter != src2.size()));

        // Insert remaining elements from src1
        if(src1_iter != src1.size()) {
            for(int i = src1_iter; i < src1.size(); i++) {
                dst[dst_iter] = src1[i];
                dst_iter++;
            }
        }

        // Insert remaining elements from src2
        if(src2_iter != src2.size()) {
            for(int i = src2_iter; i < src2.size(); i++) {
                dst[dst_iter] = src2[i];
                dst_iter++;
            }
        }

        if(num_duplicates > 0) {
            dst.resize(dst.size() - num_duplicates);
        }
    }

    // Returns any matches
    std::vector<std::pair<T, T>> insert_internal(const std::vector<T>& values)
    {
        _f = fopen(_full_path.c_str(), "rb+");
        
        std::vector<std::pair<T, T>> duplicates;

        // Sort the values
        auto sorted = values;
        std::sort(sorted.begin(), sorted.end(), GigaSet<T, CompareFunc, HashFunc>::compare_less);

        // Divide into buckets
        std::vector<T> local_buckets[BUCKETS];

        for(auto v : sorted) {
            uint32_t bucket = HashFunc(v, BUCKETS);
            local_buckets[bucket].push_back(v);
        }

        // Find the largest 
        int x1 = 0;
        int x2 = 0;
        for(int b = 0; b < BUCKETS; b++) {
            x1 = std::max(x1, (int)local_buckets[b].size());
            x2 = std::max(x2, _counts[b]);
        }

        if(x1 + x2 >= _bucket_size) {

            // Expand to next power of 2 elements
            int new_size = _bucket_size;
            while(new_size <= x1 + x2) {
                new_size *= 2;
            }
            expand_buckets(new_size);
        }

        for(int b = 0; b < BUCKETS; b++) {
            // Read from file
            std::vector<T> file_bucket(_counts[b]);
            read_bucket(file_bucket, b);
            // Merge the two lists
        
            std::vector<T> new_bucket(file_bucket.size() + local_buckets[b].size());

            merge(new_bucket, file_bucket, local_buckets[b], duplicates);

            write_bucket(b, new_bucket);
        }

        fclose(_f);

        return duplicates;
    }

    void fcopy(FILE* f, size_t src, size_t dst, size_t count)
    {
        int offset = 0;
        while(count > 0) {
            char buf[16384];

            int to_copy = count >= sizeof(buf) ? sizeof(buf) : count;
            assert(fseek(f, BASE_OFFSET + src + offset, SEEK_SET) == 0);
            assert(fread(buf, 1, to_copy, f) == to_copy);
            assert(fseek(f, BASE_OFFSET + dst + offset, SEEK_SET) == 0);
            assert(fwrite(buf, 1, to_copy, f) == to_copy);
            offset += to_copy;
            count -= to_copy;
        }
    }

    // Read entire bucket
    void read_bucket(std::vector<T>& values, int bucket)
    {
        if(_counts[bucket] > 0) {

            assert(fseek(_f, BASE_OFFSET + bucket * _bucket_size * sizeof(T), SEEK_SET) == 0);

            assert(fread(values.data(), sizeof(T), values.size(), _f) == values.size());
        }
    }

    // Write to bucket
    void write_bucket(int bucket, const std::vector<T>& data)
    {
        _counts[bucket] = data.size();

        assert(fseek(_f, bucket * sizeof(_counts[0]), SEEK_SET) == 0);
        assert(fwrite(&_counts[bucket], sizeof(_counts[0]), 1, _f) == 1);
        assert(fseek(_f, BASE_OFFSET + bucket * _bucket_size * sizeof(T), SEEK_SET) == 0);
        assert(fwrite(data.data(), sizeof(T), data.size(), _f) == data.size());
    }

    void write_counts()
    {
        assert(fseek(_f, 0, SEEK_SET) == 0);
        assert(fwrite(_counts.data(), sizeof(_counts[0]), _counts.size(), _f) == _counts.size());
    }

    void expand_buckets(int new_size)
    {
        for(int i = BUCKETS - 1; i >= 1; i--) {
            // Move bucket i from i * bucket_size to i * new_bucket_size;

            size_t src = i * _bucket_size * sizeof(T);
            size_t dst = i * new_size * sizeof(T);
            int count = _counts[i] * sizeof(T);

            fcopy(_f, src, dst, count);
        }

        _bucket_size = new_size;
    }

public:

    GigaSet(const std::string& file_path)
    {
        _full_path = file_path;

        init();
    }

    ~GigaSet()
    {
    }

    std::vector<std::pair<T, T>> insert(const std::vector<T>& values)
    {
        return insert_internal(values);
    }

    size_t size()
    {
        size_t count = 0;

        for(int b = 0; b < BUCKETS; b++) {
            count += _counts[b];
        }

        return count;
    }

};

#endif