# vecSZ

vecSZ is an implementation of the popular [SZ lossy compressor](https://github.com/szcompressor/SZ) optimized for efficient CPU performance using vectorization.

# Use
## Compilation

To compile vecSZ with default optimizations (shown below), just type `make`.

### Default Optimizations
All options are shown with their default value.
```
Enable vector support   -  vector_support=AVX2
Enable -O3 Optimization -  optimize=True
Enable OpenMP           -  omp=False
Use -g for debugging    -  debug=False
```

To enable an alternate level of vector support, use `vector_support=<None|AVX|AVX2|AVX512>`.

For example, enable AVX-512 and OpenMP by typing `make vector_support=AVX512 omp=True`.

## Basic Use

Type `vecsz` or `vecsz -h` for instant instructions. Also, a basic use vecSZ is given below.

```bash
./bin/vecsz -t f32 -m abs -e 1.0e-4.0 -i ./data/ex-cesm-CLDHGH -l 3600,1800 -z
            ------ ------ ----------- ------------------------ ------------  |
             dtype  mode  error bound      input datum file    low-to-high  zip

./bin/vecsz -i ./data/ex-cesm-CLDHGH.sz -x
            ---------------------------  |
                      sz archive        unzip
```
In this case, 1800-by-3600 (y-x order) CESM-ATM is listed in the demo set. Use `-D cesm` to specify the preset for demonstration. Type `vecsz` or `vecsz -h` to look up the presets.

```bash
vecsz -m abs -e 1e-4 -i ./data/ex-cesm-CLDHGH -D cesm -z
```
The following **essential** arguments are required,

- WORKFLOW: `-z` to compress; `-x` to decompress.
- CONFIG: `-m` to specify error control mode from `abs` (absolute) and `r2r` (relative to value range)
- CONFIG: `-e` to specify error bound
- CONFIG: `-t` to specify datatype from `f32` (single-precision float) and `f64` (double-precision float)
- INPUT: `-i` to specify input file
- INPUT: `-D` to specify demo dataset name or `-l` to dimensions
- OUTPUT: `-o` or `--output` to specify output path for both compression and decompression
- LOG: `-V` or `--verbose` to print CPU information

# Limitations
- Current support for only x86 architectures (AVX/SSE support)
- Support for single-precision float data only (f32)
