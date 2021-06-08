# vecSZ

vecSZ is an implementation of the popular [SZ lossy compressor](https://github.com/szcompressor/SZ) optimized for efficient CPU performance using vectorization.

# Use
## basic use

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
- INPUT: `-i` to specify input file
- INPUT: `-D` to specify demo dataset name or `-l` to dimensions
- OUTPUT: `-o` or `--output` to specify output path for both compression and decompression
- LOG: `-V` or `--verbose` to print host and device information

# Limitations
- Current support for only x86 architectures (AVX/SSE support)
- Support for single-precision float data only (f32)
- Block size must be a multiple of the vector register used (e.g. 8, 16, 32, 64)
