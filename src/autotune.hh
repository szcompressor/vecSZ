#ifndef PSZ_TUNING_HH
#define PSZ_TUNING_HH

#include <iostream>
#include "psz_dualquant_opt.hh"
#include "argument_parser/argparse.hh"
#include "types.hh"
#include "dimensions.hh"

/* Used for timing */
using namespace std;
using namespace std::chrono;

namespace PdQ  = pSZ::PredictionDualQuantization;

high_resolution_clock::time_point timept = high_resolution_clock::now();
#define TIMER duration_cast<duration<double>>(high_resolution_clock::now() - timept).count()

/* TODO: - TUNE BLOCK SIZE          */
/*       - TUNE VECTORIZATION LEVEL */
/*       - TUNE PREFETCHING         */

template <typename T, typename Q>
double run_sample_blocks(argparse* ap, bool fine_massive, bool blocked, T* data, T* outlier, Q* code, size_t const* const dims, double const* const ebs_L4, T* pred_err, T* comp_err)
{
    auto num_iterations = ap->num_iterations;
    auto sample_pct     = ap->sample_percentage;
    auto blksz          = ap->block_size;
    auto vecsz          = ap->vector_length;

    //cout << "Testing Block Size = " << blksz << ", Vector Registers = " << vecsz << endl;
    auto dims_L16 = InitializeDims(ap);

    size_t nblocks  = (dims_L16[nBLK0] * dims_L16[nBLK1] * dims_L16[nBLK2] * dims_L16[nBLK3]);
    size_t nsamples = (sample_pct * 0.01) * nblocks + 1;
    //cout << "Number of Samples: " << nsamples << endl;

    int iterations = num_iterations;
    for (int i = 0; i < iterations + 1; i++) {
	srand(1);
    if (i == 1) timept = high_resolution_clock::now(); //begin timing

    if (dims_L16[nDIM] == 1) {
        if (blocked) {
#pragma omp parallel for
            for (size_t n = 0; n < nsamples; n++) {
                size_t b0 = rand() % dims_L16[nBLK0];
                if (fine_massive)
#ifdef REPBLK
                    PdQ::c_lorenzo_1d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, 0, blksz, vecsz);
#else
                    PdQ::c_lorenzo_1d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, b0, blksz, vecsz);
#endif
            }
        }
    } else if (dims_L16[nDIM] == 2) {
        if (blocked) {
#pragma omp parallel for
            for (size_t n = 0; n < nsamples; n++) {
                size_t b0 = rand() % dims_L16[nBLK0];
                size_t b1 = rand() % dims_L16[nBLK1];
                if (fine_massive)
#ifdef REPBLK
                    PdQ::c_lorenzo_2d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, 0, 0, blksz, vecsz);
#else
                    PdQ::c_lorenzo_2d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, b0, b1, blksz, vecsz);
#endif
            }
        } 
    } else if (dims_L16[nDIM] == 3) {
        if (blocked) {
#pragma omp parallel for
            for (size_t n = 0; n < nsamples; n++) {
                size_t b0 = rand() % dims_L16[nBLK0];
                size_t b1 = rand() % dims_L16[nBLK1];
                size_t b2 = rand() % dims_L16[nBLK2];
                if (fine_massive)
#ifdef REPBLK
                    PdQ::c_lorenzo_3d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, 0, 0, 0, blksz, vecsz);
#else
                    PdQ::c_lorenzo_3d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, b0, b1, b2, blksz, vecsz);
#endif
            }
        } 
    }
    }//end iterations
    double sample_time = TIMER / (iterations); //end timing

    //printf("\tAverage Sample Duration: %0.6f\n", sample_time);

    return sample_time;
}

template <typename T, typename Q>
int autotune_block_sizes(argparse* ap, double* timing, bool fine_massive, bool blocked, T* data, T* outlier, Q* bcode, size_t const* const dims, double const* const ebs_L4)
{
	const int NSIZES = 4;
	const int sizes[NSIZES] = {8, 16, 32, 64};
	int blksz;
	double newTime;

	*timing = run_sample_blocks<T, Q>(ap, fine_massive, blocked, data, outlier, bcode, dims, ebs_L4);
	for (int i = 1; i < NSIZES; i++)
	{
		newTime = run_sample_blocks<T, Q>(ap, fine_massive, blocked, data, outlier, bcode, dims, ebs_L4);
		blksz = (newTime < *timing) ? sizes[i] : blksz;
		*timing = (newTime < *timing) ? newTime : *timing;
	}

	return blksz;
}

template <typename T, typename Q>
int autotune_vector_len(argparse* ap, int* blksz, double* timing, bool fine_massive, bool blocked, T* data, T* outlier, Q* bcode, size_t const* const dims, double const* const ebs_L4)
{
	const int NSIZES = 2;
	const int sizes[NSIZES] = {256, 512};
	int vecsz, blksz_256, blksz_512;
	double time_512, time_256;

	blksz_256 = autotune_block_sizes<T,Q>(ap, &time_256, fine_massive, blocked, data, outlier, bcode, dims, ebs_L4);
#ifdef AVX512
	blksz_512 = autotune_block_sizes<T,Q>(ap, &time_512, fine_massive, blocked, data, outlier, bcode, dims, ebs_L4);
	
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
     
    cout << "Best Performance: block_size = " << *blksz << ", vector_length = " << vecsz << ", sample_pct = " << ap->sample_percentage << ", num_iters = " << ap->num_iterations << endl;  
	
	if (*blksz > 64 || *blksz < 8) *blksz = 64;
	return vecsz;
}


#endif
