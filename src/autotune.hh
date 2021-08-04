#ifndef PSZ_TUNING_HH
#define PSZ_TUNING_HH

#include <iostream>
#include "dualquant.hh"
#include "argument_parser/argparse.hh"
#include "types.hh"
#include "dimensions.hh"
#include "utils/timer.hh"

/* Used for timing */
using namespace std;
using namespace std::chrono;

namespace pq  = vecsz::predictor_quantizer;

/* TODO: - TUNE BLOCK SIZE          */
/*       - TUNE VECTORIZATION LEVEL */
/*       - TUNE PREFETCHING         */

template <typename T, typename Q>
double run_sample_blocks(argparse* ap, T* data, T* outlier, Q* code, size_t const* const dims, double const* const ebs_L4)
{
    hires::time_point start;
    auto num_iterations = ap->num_iterations;
    auto sample_pct     = ap->sample_percentage;
    auto blksz          = ap->block_size;
    auto vecsz          = ap->vector_length;
    auto dims_L16       = InitializeDims(ap);

    size_t nblocks  = (dims_L16[nBLK0] * dims_L16[nBLK1] * dims_L16[nBLK2] * dims_L16[nBLK3]);
    size_t nsamples = (sample_pct / 100) * nblocks + 1;

    int iterations = num_iterations;
    for (int i = 0; i < iterations + 1; i++)
    {
        srand(1);
        if (i == 1) start = hires::now(); //begin timing

        if (dims_L16[nDIM] == 1)
        {
            #pragma omp parallel for
            for (size_t n = 0; n < nsamples; n++)
            {
                size_t b0 = rand() % dims_L16[nBLK0];
                pq::c_lorenzo_1d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, b0, blksz, vecsz);
            }
        } else if (dims_L16[nDIM] == 2)
        {
            #pragma omp parallel for
            for (size_t n = 0; n < nsamples; n++)
            {
                size_t b0 = rand() % dims_L16[nBLK0];
                size_t b1 = rand() % dims_L16[nBLK1];
                pq::c_lorenzo_2d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, b0, b1, blksz, vecsz, ap->szwf, ap->pad_constant, ap->pad_type);
            }
        }
        else if (dims_L16[nDIM] == 3)
        {
            #pragma omp parallel for
            for (size_t n = 0; n < nsamples; n++)
            {
                size_t b0 = rand() % dims_L16[nBLK0];
                size_t b1 = rand() % dims_L16[nBLK1];
                size_t b2 = rand() % dims_L16[nBLK2];
                pq::c_lorenzo_3d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, b0, b1, b2, blksz, vecsz, ap->szwf, ap->pad_constant, ap->pad_type);
            }
        }
    }//end iterations

    auto end = hires::now();
    return static_cast<duration_t>(end - start).count() / (iterations); //end timing
}

template <typename T, typename Q>
int autotune_block_sizes(argparse* ap, double* timing, T* data, T* outlier, Q* bcode, size_t const* const dims, double const* const ebs_L4)
{
	const int NSIZES = 4;
	const int sizes[NSIZES] = {8, 16, 32, 64};
	int blksz;
	double newTime;

	*timing = run_sample_blocks<T, Q>(ap, data, outlier, bcode, dims, ebs_L4);
	for (int i = 1; i < NSIZES; i++)
	{
		newTime = run_sample_blocks<T, Q>(ap, data, outlier, bcode, dims, ebs_L4);
		blksz = (newTime < *timing) ? sizes[i] : blksz;
		*timing = (newTime < *timing) ? newTime : *timing;
	}

	return blksz;
}

template <typename T, typename Q>
int autotune_vector_len(argparse* ap, int* blksz, double* timing, T* data, T* outlier, Q* bcode, size_t const* const dims, double const* const ebs_L4)
{
	const int NSIZES = 2;
	const int sizes[NSIZES] = {256, 512};
	int vecsz, blksz_256, blksz_512;
	double time_512, time_256;

	blksz_256 = autotune_block_sizes<T,Q>(ap, &time_256, data, outlier, bcode, dims, ebs_L4);
#ifdef AVX512
	blksz_512 = autotune_block_sizes<T,Q>(ap, &time_512, data, outlier, bcode, dims, ebs_L4);

	if (time_512 < time_256) {
		*blksz = blksz_512;
		vecsz = 512;
	}
	else {
		*blksz = blksz_256;
		vecsz = 256;
	}
#else
	*blksz = blksz_256;
	vecsz = 256;
#endif

    LogAll(log_cfg, "sample percentage:", ap->sample_percentage,"number repetitions:",ap->num_iterations);
    LogAll(log_cfg, "best performance --", "block size:",*blksz,"vector length:",vecsz);

	if (*blksz > 64 || *blksz < 8) *blksz = 64;
	return vecsz;
}


#endif
