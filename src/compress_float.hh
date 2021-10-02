#ifndef COMPRESS_FLOAT_HH
#define COMPRESS_FLOAT_HH

#include <cstddef>
#include <cstdlib>

#include "constants.hh"
#include "huffman.hh"
#include "lossless.hh"
#include "argument_parser/argparse.hh"

/* namespace DV = DesignVerification; */

template <typename T>
T max(T* array, size_t len)
{
    T max = array[0];
    for (size_t i = 0; i < len; i++)
    {
        max = (array[i] > max) ? array[i] : max;
    }

    return max;
}

template <typename T>
T min(T* array, size_t len)
{
    T min = array[0];
    for (size_t i = 0; i < len; i++)
    {
        min = (array[i] < min) ? array[i] : min;
    }

    return min;
}

template <typename T> T median(T*array, size_t len) { return (max(array, len) - min(array, len)) / 2; } namespace vecsz {

template <typename T>
unsigned char* compress_float_arr(SZWorkflow szwf, T* array, unsigned int length, unsigned int* compressed_length)
{
    int* iarray = new int[length];
    for (unsigned int i = 0; i < length; i++) iarray[i] = static_cast<int>(array[i]);
    unsigned char* new_array = new unsigned char[length * sizeof(int)];
    memcpy(new_array, iarray, sizeof(int) * length);
    delete[] iarray;

	unsigned char* lossless_out = nullptr;
    *compressed_length = 0;

    if (szwf.lossless_gzip) *compressed_length = vecsz::lossless::sz_lossless_compress(GZIP_COMPRESSOR, 3, new_array, length * sizeof(int), &lossless_out);
    else *compressed_length = vecsz::lossless::sz_lossless_compress(ZSTD_COMPRESSOR, 3, new_array, length * sizeof(int), &lossless_out);

    delete[] new_array;
    return lossless_out;
}

template <typename T>
T* decompress_float_arr(unsigned char* array, size_t length, size_t target_size)
{

    unsigned char *data = NULL;
    int lossless_pass = vecsz::lossless::is_lossless_compressed_data(array, length);
    vecsz::lossless::sz_lossless_decompress(lossless_pass, array, length, &data, target_size * sizeof(int));

    int* new_array = new int[target_size];
    memcpy(new_array, data, sizeof(int) * target_size);

    T* output_array = new T[target_size];
    for (unsigned int i = 0; i < target_size; i++) output_array[i] = static_cast<T>(new_array[i]);

    return output_array;
}

} // vecsz

#endif /* COMPRESS_FLOAT_HH */
