/**
 * @file document.hh
 * @author Jiannan Tian
 * @brief
 * @version 0.1.1
 * @date 2020-09-22
 *
 * @copyright (C) 2020 by Washington State University, Argonne National Laboratory
 * See LICENSE in top-level directory
 *
 */

#ifndef ARGUMENT_PARSER_DOCUMENT_HH
#define ARGUMENT_PARSER_DOCUMENT_HH

#include <regex>
#include <string>

using std::regex;
using std::string;

const string fmt_b("\e[1m");
const string fmt_0("\e[0m");

const regex  bful("@(.*?)@");
const string bful_text("\e[1m\e[4m$1\e[0m");
const regex  bf("\\*(.*?)\\*");
const string bf_text("\e[1m$1\e[0m");
const regex  ul(R"(_((\w|-|\d|\.)+?)_)");
const string ul_text("\e[4m$1\e[0m");
const regex  red(R"(\^\^(.*?)\^\^)");
const string red_text("\e[31m$1\e[0m");

string  //
Format(const std::string& s)
{
    auto a = std::regex_replace(s, bful, bful_text);
    auto b = std::regex_replace(a, bf, bf_text);
    auto c = std::regex_replace(b, ul, ul_text);
    auto d = std::regex_replace(c, red, red_text);
    return d;
}

static const char vecsz_short_doc[] =
    // "vecsz, version [placeholder]\n"
    "\n"
    "usage: vecsz [-zxrhV] [-i file] [-t dtype] [-m mode] [-e eb] [-l x,y,z] ...\n"
    "\n"
    " -z : zip/compress\n"
    " -x : unzip/decompress\n"
    " -r : dryrun\n"
    " -h : print full-length help document\n"
    " -V : be verbose\n"
    "\n"
    " -i file        : path to input datum\n"
    " -t dtype       : f32\n" //TODO: f[32|64] or fp[4|8], i[8|16|32|64] or int[1|2|3|4]\n"
    " -m mode        : compression mode; abs, r2r\n"
    " -e eb          : error bound; default 1e-4\n"
    " -l [x[,y[,z]]] : Specify (1|2|3)D datum/field sizes, with dimensions from low to high.\n"
    " -D name        : use demo dataset, skip interpretation\n"
    "                    (1D) hacc  hacc1b  (2D) cesm  exafel\n"
    "                    (3D) hurricane  nyx-s  nyx-m  qmc  qmcpre  aramco  parihaka\n"
    " -d num         : specify codebook size\n"
    " -s vector_len  : specify max size vector register to use in bits or (AVX, AVX2, AVX512, etc.)\n"
    " -b block_size  : specify block size to use for chunking, default 8\n"
    " -S to_skip     : specify portions of the vecSZ workflow to skip (huffman, write, verify)\n"
    " -L [gzip|zstd] : perform lossless pass with either gzip or zstd\n"
    "\n"
    "example:\n"
    "  zip 1: ./bin/vecsz -t f32 -m r2r -e 1e-4 -i ./data/ex-cesm-CLDHGH -D cesm -z\n"
    "  zip 2: ./bin/vecsz -t f32 -m r2r -e 1e-4 -i ./data/ex-cesm-CLDHGH -l 3600,1800 -z\n"
    "  unzip: ./bin/vecsz -i ./data/ex-cesm-CLDHGH.sz -x\n"
    "\n"
    "\"vecsz -h\" for details.\n";

static const char vecsz_full_doc[] =
    "*NAME*\n"
    "        vecSZ: Efficient Error-Bounded Lossy Compressor for Scientific Data\n"
    "        Lowercased \"*vecsz*\" is the command."
    "\n"
    "*SYNOPSIS*\n"
    "        The basic use is listed below,\n"
    "        *vecsz* *-t* f32 *-m* r2r *-e* 1.0e-4.0 *-i* ./data/ex-cesm-CLDHGH *-l* 3600,1800 *-z* \n"
    //       vecsz -t f32 -m r2r -e 1.0e-4.0 -i ./data/ex-cesm-CLDHGH -l 3600,1800 -z \n
    "             ------ ------ ----------- ------------------------ ------------  | \n"
    "              dtype  mode  error bound      input datum file    low-to-high  zip \n"
    "\n"
    "        *vecsz* *-i* ./data/ex-cesm-CLDHGH.sz *-x* \n"
    //       vecsz -i ./data/ex-cesm-CLDHGH.sz -x \n"
    "             ---------------------------  | \n"
    "                       sz archive        unzip \n"
    "\n"
    "        *vecsz* *-t* f32|64 *-m* [eb mode] *-e* [eb] *-i* [datum file] *-D* [demo dataset] *-z*\n"
    "        *vecsz* *-t* f32|64 *-m* [eb mode] *-e* [eb] *-i* [datum file] *-l* [nx[,ny[,nz]]] *-z*\n"
    "        *vecsz* *-i* [datum basename].sz *-x*\n"
    "\n"
    "*OPTIONS*\n"
    "    *Mandatory* (zip and dryrun)\n"
    "        *-z* or *--compress* or *--*@z@*ip*\n"
    "        *-r* or *--dry-*@r@*un*\n"
    "                Simulate compress/decompress workflow without writing output.\n"
    "                Timing for complete workflow and data quality report is produced.\n"
    "\n"
    "        *-i* or *--*@i@*nput* [datum file]\n"
    "\n"
    "        *-l* [x[,y[,z]]]\n"
    "                Specify (1|2|3)D datum/field sizes, with dimensions from low to high.\n"
    "\n"
    "        *-m* or *--*@m@*ode* <abs|r2r>\n"
    "                Specify error-controlling mode. Supported modes include:\n"
    "                _abs_: absolute mode, eb = input eb\n"
    "                _r2r_: relative-to-value-range mode, eb = input eb x value range\n"
    "\n"
    "        *-e* or *--eb* or *--error-bound* [num]\n"
    "                Specify error bound. e.g., _1.23_, _1e-4_, _1.23e-4.56_\n"
    "\n"
    "        *-d* or *--dict-size* [256|512|1024|...]\n"
    "                Specify dictionary size/quantization bin number.\n"
    "                Should be a power-of-2.\n"
    "\n"
    "        *-b* or *--block-size* [num]\n"
    "                Specify block size used for chunking data.\n"
    "\n"
    "        *-s* or *--vector* [nbits|AVX|AVX2|AVX512]\n"
    "                Specify max vector register sizes in bits or by vector extension name.\n"
    "\n"
    "        *-S* or *--*@s@*kip* [write|huffman|verify]\n"
    "                Skip one or more specified steps in the vecSZ workflow.\n"
    "                For skipping multiple steps use a comma to separate each one.\n"
    "\n"
    "        *-L* or *--*@l@*ossless* [gzip|zstd]\n"
    "                Perform a lossless pass using the specified lossless compressor,\n"
    "                GZIP or ZSTD.\n"
    "\n"
    "        *-a* or *--*@a@*utotune*\n"
    "                Perform autotuning for block size and vector length before\n"
    "                beginning compress/decompress.\n"
    "\n"
    "        *--num-iter* [num]\n"
    "                For autotuning, choose the number of iterations to run tuning\n" 
    "                before selecting an optimal configuration. Default 10.\n"
    "\n"
    "        *--sample* [num]\n"
    "                For autotuning, choose the percentage of the total dataset (0-100)\n" 
    "                to test when tuning. Default: 1\%.\n"
    "\n"
    "    *Mandatory* (unzip)\n"
    "        *-x* or *--e*@x@*tract* or *--decompress* or *--unzip*\n"
    "\n"
    "        *-i* or *--*@i@*nput* [corresponding datum basename (w/o extension)]\n"
    "\n"
    "    *Additional I/O*\n"
    "        *--input* /path/to/origin-datum\n"
    "                For verification & get data quality evaluation.\n"
    "        *--output*  /path/to\n"
    "                Specify alternative output path.\n"
    "\n"
/*    "    *Modules*\n"
    "        *-S* or *--e*@x@*clude* or *--*@s@*kip* _module-1_,_module-2_,...,_module-n_,\n"
    "                Disable functionality modules. Supported module(s) include:\n"
    "                _huffman_  Huffman codec after prediction+quantization (p+q) and before reversed p+q.\n"
    "                _write.x_  Skip write decompression data.\n"
    "\n"
    "        *-p* or *--pre* _method-1_,_method-2_,...,_method-n_\n"
    "                Enable preprocessing. Supported preprocessing method(s) include:\n"
    "                _binning_  Downsampling datum by 2x2 to 1.\n"
    "\n"
*/    "    *Demonstration*\n"
    "        *-h* or *--help*\n"
    "                Get help documentation.\n"
    "\n"
    "        *-V* or *--verbose*\n"
    "                Print host and device information for diagnostics.\n"
    "\n"
 /*   "        *-M* or *--meta*\n"
    "                Get archive metadata. (TODO)\n"
    "\n"
 */   "        *-D* or *--demo* [demo-dataset]\n"
    "                Use demo dataset, will omit given dimension(s). Supported datasets include:\n"
    "                1D: _hacc_  _hacc1b_    2D: _cesm_  _exafel_\n"
    "                3D: _hurricane_  _nyx-s_  _nyx-m_  _qmc_  _qmcpre_  _aramco_  _parihaka_\n"
    "\n"
  /*  "    *Internal* (will be automated with configuration when going public)\n"
    "        *-Q* or *--*@q@*uant-byte* <1|2>\n"
    "                Specify quantization code representation.\n"
    "                Options _1_, _2_ are for *1-* and *2-*byte, respectively. (default: 2)\n"
    "                ^^Manually specifying this may not result in optimal memory footprint.^^\n"
    "\n"
    "        *-H* or *--*@h@*uff-byte* <4|8>\n"
    "                Specify Huffman codeword representation.\n"
    "                Options _4_, _8_ are for *4-* and *8-*byte, respectively. (default: 4)\n"
    "                ^^Manually specifying this may not result in optimal memory footprint.^^\n"
    "\n"
    "        *-C* or *--huff-*@c@*hunk* [256|512|1024|...]\n"
    "                Manually specify chunk size for Huffman codec, overriding autotuning.\n"
    "                Should be a power-of-2 that is sufficiently large.\n"
    "                ^^This affects Huffman decoding performance significantly.^^\n"
    "\n"
*/    "*EXAMPLES*\n"
    "    *Demo Datasets*\n"
    "        *CESM* example:\n"
    "        ./bin/vecsz -t f32 -m r2r -e 1e-4 -i ./data/ex-cesm-CLDHGH -D cesm -z\n"
    "        ./bin/vecsz -t f32 -m r2r -e 1e-4 -i ./data/ex-cesm-CLDHGH -D cesm -r\n"
    "        ./bin/vecsz -i ./data/ex-cesm-CLDHGH.sz -x\n"
    "\n"
/*    "        *CESM* example with specified output path:\n"
    "        mkdir data2 data3\n"
    "            # zip, output to `data2`\n"
    "        ./bin/vecsz -t f32 -m r2r -e 1e-4 -i ./data/ex-cesm-CLDHGH -D cesm -z --output data2\n"
    "            # unzip, in situ\n"
    "        ./bin/vecsz -i ./data2/ex-cesm-CLDHGH.sz -x && ls data2\n"
    "            # unzip, output to `data3`\n"
    "        ./bin/vecsz -i ./data2/ex-cesm-CLDHGH.sz -x --output data3 && ls data3\n"
    "            # unzip, output to `data3`, compare to the original datum\n"
    "        ./bin/vecsz -i ./data2/ex-cesm-CLDHGH.sz -x --output data3 --origin ./data/ex-cesm-CLDHGH && ls "
    "data3\n"
    "        ## Please create directory by hand before specifying (considering access permission control).\n"
    "\n"
 */   "        *Hurricane Isabel* example:\n"
    "        ./bin/vecsz -t f32 -m r2r -e 1e-4 -i ./data/ex-hurr-CLOUDf48 -D hurricane -z\n"
    "        ./bin/vecsz -t f32 -m r2r -e 1e-4 -i ./data/ex-hurr-CLOUDf48 -D hurricane -r\n"
    "        ./bin/vecsz -i ./data/ex-hurr-CLOUDf48.sz -x\n"
    "\n";

static const char doc_dim_order[] =
    "\n"
    "  Input dimension follows low-to-high (e.g., x-y-z) order.\n"
    "  Taking 2D CESM-ATM as an example, \n"
    "\n"
    "  |<------------------------- x 3600 --------------------------->|    \n"
    "  +--------------------------------------------------------------+  - \n"
    "  |                                                              |  ^ \n"
    "  |                                                              |  | \n"
    "  |              CESM-ATM:    1800x3600 (y-x order)              |  | \n"
    "  |              datum name:  <field>_1800_3600                  |  y \n"
    "  |                                                              | 1800 \n"
    "  |              input:       -l 3600,1800                       |  | \n"
    "  |              input order: -l [x,y]                           |  | \n"
    "  |                                                              |  | \n"
    "  |                                                              |  v \n"
    "  +--------------------------------------------------------------+  - \n"
    "\n"
    "  Taking 3D Hurricane as another example, whose dimensions are\n"
    "  100x500x500, the input is \"-l 500,500,100\".\n";

#endif
