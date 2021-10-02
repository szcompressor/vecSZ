#ifndef LOSSLESS_HH
#define LOSSLESS_HH
/**
 * @file lossless.hh
 * @author Griffin Dube
 * @brief gzip and zstd compressor interface: adapted from SZ
 *             callZlib.c and utility.c by Sheng Di
 * @version 0.1
 * @date 2021-07-01
 *
 * @copyright Copyright (c) 2021
 *
 */

#include <iostream>
#include <cstdlib>
#include <zlib.h>

#include "constants.hh"
#include "utils/format.hh"
#include "zstd/zstd.h"
#include "zlib/zlib.h"

using std::cerr;
using std::cout;
using std::endl;

#if MAX_MEM_LEVEL >= 8
#define DEF_MEM_LEVEL 8
#else
#define DEF_MEM_LEVEL MAX_MEM_LEVEL
#endif

#define SZ_ZLIB_BUFFER_SIZE 65536
#define GZIP_COMPRESSOR 0
#define ZSTD_COMPRESSOR 1

namespace vecsz {

    namespace lossless {

        int isZlibFormat(unsigned char magic1, unsigned char magic2)
        {
            if(magic1==104&&magic2==5) //DC+BS
                return 1;
            if(magic1==104&&magic2==129) //DC+DC
                return 1;
            if(magic1==104&&magic2==222) //DC+BC
                return 1;
            if(magic1==120&&magic2==1) //BC+BS
                return 1;
            if(magic1==120&&magic2==94) //BC+?
                return 1;
            if(magic1==120&&magic2==156) //BC+DC
                return 1;
            if(magic1==120&&magic2==218) //BC+BS
                return 1;
            return 0;
        }

        int is_lossless_compressed_data(unsigned char* compressedBytes, size_t cmpSize)
        {
        #if ZSTD_VERSION_NUMBER >= 10300
            unsigned long long frameContentSize = ZSTD_getFrameContentSize(compressedBytes, cmpSize);
            if(frameContentSize != ZSTD_CONTENTSIZE_ERROR)
                return ZSTD_COMPRESSOR;
        #else
            unsigned long long frameContentSize = ZSTD_getDecompressedSize(compressedBytes, cmpSize);
            if(frameContentSize != 0)
                return ZSTD_COMPRESSOR;
        #endif
            int flag = isZlibFormat(compressedBytes[0], compressedBytes[1]);
            if(flag)
                return GZIP_COMPRESSOR;

            return -1; //fast mode (without GZIP or ZSTD)
        }

        /* zlib_compress5() from SZ */
        unsigned long zlib_compress(unsigned char* data, unsigned long dataLength, unsigned char** compressBytes, int level)
        {
            int ret, flush;
            unsigned have;
            z_stream strm;
            unsigned char* in = data;

            /* allocate deflate state */
            strm.zalloc = Z_NULL;
            strm.zfree = Z_NULL;
            strm.opaque = Z_NULL;
            ret = deflateInit(&strm, level);
            //int windowBits = 15;
            //ret = deflateInit2(&strm, level, Z_DEFLATED, windowBits, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);//Z_FIXED); //Z_DEFAULT_STRATEGY

            if (ret != Z_OK)
                return ret;

            size_t p_size = 0, av_in = 0;
            uLong estCmpLen = deflateBound(&strm, dataLength);
            *compressBytes = (unsigned char*)malloc(sizeof(unsigned char)*estCmpLen);
            unsigned char* out = *compressBytes;

            /* compress until end of file */
            do {
                p_size += SZ_ZLIB_BUFFER_SIZE;
                if(p_size>=dataLength) { av_in = dataLength - (p_size - SZ_ZLIB_BUFFER_SIZE); flush = Z_FINISH;
                }
                else
                {
                    av_in = SZ_ZLIB_BUFFER_SIZE;
                    flush = Z_NO_FLUSH;
                }
                strm.avail_in = av_in;
                strm.next_in = in;

                /* run deflate() on input until output buffer not full, finish
                compression if all of source has been read in */
                do {
                    strm.avail_out = SZ_ZLIB_BUFFER_SIZE;
                    strm.next_out = out;
                    ret = deflate(&strm, flush);    /* no bad return value */

                    have = SZ_ZLIB_BUFFER_SIZE - strm.avail_out;
                    out += have;
                } while (strm.avail_out == 0);

                in+=av_in;

                /* done when last data in file processed */
            } while (flush != Z_FINISH);

            /* clean up and return */
            (void)deflateEnd(&strm);

            return strm.total_out;
        }

        /* zlib_uncompress5() from SZ */
        unsigned long zlib_uncompress(unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize)
        {
            int err;
            z_stream d_stream = {0}; /* decompression stream */

            *oriData = (unsigned char*)malloc(sizeof(unsigned char)*targetOriSize);

            d_stream.zalloc = (alloc_func)0;
            d_stream.zfree = (free_func)0;
            d_stream.opaque = (voidpf)0;

            d_stream.next_in  = compressBytes;
            d_stream.avail_in = 0;
            d_stream.next_out = *oriData;

            err = inflateInit(&d_stream);
            if (err != Z_OK && err != Z_STREAM_END)
            {
                cerr << log_err << "inflateInit error: " << err << endl;
            }

            while (d_stream.total_out < targetOriSize && d_stream.total_in < cmpSize) {
                d_stream.avail_in = d_stream.avail_out = SZ_ZLIB_BUFFER_SIZE; /* force small buffers */
                //err = inflate(&d_stream, Z_NO_FLUSH);
                err = inflate(&d_stream, Z_SYNC_FLUSH);
                if (err == Z_STREAM_END) break;
                if (err != Z_OK && err != Z_STREAM_END)
                {
                    cerr << log_err << "inflate error: " << err << endl;
                }
            }

            err = inflateEnd(&d_stream);

            if (err != Z_OK && err != Z_STREAM_END)
            {
                cerr << log_err << "inflateEnd error: " << err << endl;
            }

            return d_stream.total_out;
        }

        unsigned long sz_lossless_compress(int losslessCompressor, int level, unsigned char* data, unsigned long dataLength, unsigned char** compressBytes)
        {
            unsigned long outSize = 0;
            size_t estimatedCompressedSize = 0;
            switch(losslessCompressor)
            {
            case GZIP_COMPRESSOR:
                outSize = zlib_compress(data, dataLength, compressBytes, level);
                break;
            case ZSTD_COMPRESSOR:
                if(dataLength < 100)
                    estimatedCompressedSize = 200;
                else
                    estimatedCompressedSize = dataLength*1.2;
                *compressBytes = (unsigned char*)malloc(estimatedCompressedSize);
                outSize = ZSTD_compress(*compressBytes, estimatedCompressedSize, data, dataLength, level); //default setting of level is 3
                break;
            default:
                printf("Error: Unrecognized lossless compressor in sz_lossless_compress()\n");
            }
            return outSize;
        }

        unsigned long sz_lossless_decompress(int losslessCompressor, unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize)
        {
            unsigned long outSize = 0;
            switch(losslessCompressor)
            {
            case GZIP_COMPRESSOR:
                outSize = zlib_uncompress(compressBytes, cmpSize, oriData, targetOriSize);
                break;
            case ZSTD_COMPRESSOR:
                *oriData = (unsigned char*)malloc(targetOriSize);
                ZSTD_decompress(*oriData, targetOriSize, compressBytes, cmpSize);
                outSize = targetOriSize;
                break;
            default:
                printf("Error: Unrecognized lossless compressor in sz_lossless_decompress()\n");
            }
            return outSize;
        }
    } // namespace lossless

} // namespace vecsz

#endif /* ifndef LOSSLESS_HH */
