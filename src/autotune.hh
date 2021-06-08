#ifndef PSZ_TUNING_HH
#define PSZ_TUNING_HH

#include <iostream>
#include "SDRB.hh"
#include "psz_14.hh"
#include "psz_14blocked.hh"
#include "psz_dualquant_opt.hh"

/* Used for timing */
using namespace std;
using namespace std::chrono;

namespace PdQ  = pSZ::PredictionDualQuantization;
namespace PQRs = pSZ::PredictionQuantizationReconstructionSingleton;
namespace PQRb = pSZ::PredictionQuantizationReconstructionBlocked;

high_resolution_clock::time_point timept = high_resolution_clock::now();
#define TIMER duration_cast<duration<double>>(high_resolution_clock::now() - timept).count()

/* TODO: - TUNE BLOCK SIZE          */
/*       - TUNE VECTORIZATION LEVEL */
/*       - TUNE PREFETCHING         */

template <typename T, typename Q, size_t B>
double run_sample_blocks(int num_iterations, float sample_pct, const int blksz, const int vecsz, bool fine_massive, bool blocked, std::string const& dataset, const size_t dict_size, T* data, T* outlier, Q* code, size_t const* const dims, double const* const ebs_L4, T* pred_err, T* comp_err)
{
    //cout << "Testing Block Size = " << blksz << ", Vector Registers = " << vecsz << endl;
    auto dims_L16 = InitializeDemoDims(dataset, dict_size, blksz);

    size_t nblocks  = (dims_L16[nBLK0] * dims_L16[nBLK1] * dims_L16[nBLK2] * dims_L16[nBLK3]);
    size_t nsamples = sample_pct * nblocks + 1;
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
                    PdQ::c_lorenzo_1d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, 0, blksz, vecsz);
#else
                    PdQ::c_lorenzo_1d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0, blksz, vecsz);
#endif
                else
                    PQRb::c_lorenzo_1d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0);
            }
        } else {
            PQRs::c_lorenzo_1d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err);
        }
    } else if (dims_L16[nDIM] == 2) {
        if (blocked) {
#pragma omp parallel for
            for (size_t n = 0; n < nsamples; n++) {
                size_t b0 = rand() % dims_L16[nBLK0];
                size_t b1 = rand() % dims_L16[nBLK1];
                if (fine_massive)
#ifdef REPBLK
                    PdQ::c_lorenzo_2d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, 0, 0, blksz, vecsz);
#else
                    PdQ::c_lorenzo_2d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0, b1, blksz, vecsz);
#endif
                else
                    PQRb::c_lorenzo_2d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0, b1);
                
            }
        } else {
            PQRs::c_lorenzo_2d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err);
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
                    PdQ::c_lorenzo_3d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, 0, 0, 0, blksz, vecsz);
#else
                    PdQ::c_lorenzo_3d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0, b1, b2, blksz, vecsz);
#endif
                else
                    PQRb::c_lorenzo_3d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0, b1, b2);
            }
        } else {
            PQRs::c_lorenzo_3d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err);
        }
    }
    }//end iterations
    double sample_time = TIMER / (iterations); //end timing

    //printf("\tAverage Sample Duration: %0.6f\n", sample_time);

    return sample_time;
}

template <typename T, typename Q, size_t B>
int autotune_block_sizes(int num_iter, float sample_pct, int vec_type, double* timing, bool fine_massive, bool blocked, std::string const& dataset, const size_t dict_size, T* data, T* outlier, Q* bcode, size_t const* const dims, double const* const ebs_L4, T* pred_err, T* comp_err)
{
	const int NSIZES = 4;
	const int sizes[NSIZES] = {8, 16, 32, 64};
	int blksz;
	double newTime;

	*timing = run_sample_blocks<T, Q, B>(num_iter, sample_pct, (const int)sizes[0], vec_type, fine_massive, blocked, dataset, dict_size, data, outlier, bcode, dims, ebs_L4, pred_err, comp_err);
	for (int i = 1; i < NSIZES; i++)
	{
		newTime = run_sample_blocks<T, Q, B>(num_iter, sample_pct, (const int)sizes[i], vec_type, fine_massive, blocked, dataset, dict_size, data, outlier, bcode, dims, ebs_L4, pred_err, comp_err);
		blksz = (newTime < *timing) ? sizes[i] : blksz;
		*timing = (newTime < *timing) ? newTime : *timing;
	}

	return blksz;
}

template <typename T, typename Q, size_t B>
int autotune_vector_len(int num_iter, float sample_pct, int* blksz, double* timing, bool fine_massive, bool blocked, std::string const& dataset, const size_t dict_size, T* data, T* outlier, Q* bcode, size_t const* const dims, double const* const ebs_L4, T* pred_err, T* comp_err)
{
	const int NSIZES = 2;
	const int sizes[NSIZES] = {256, 512};
	int vecsz, blksz_256, blksz_512;
	double time_512, time_256;

	blksz_256 = autotune_block_sizes<T,Q,B>(num_iter, sample_pct, sizes[0], &time_256, fine_massive, blocked, dataset, dict_size, data, outlier, bcode, dims, ebs_L4, pred_err, comp_err);
#ifdef AVX512
	blksz_512 = autotune_block_sizes<T,Q,B>(num_iter, sample_pct, sizes[1], &time_512, fine_massive, blocked, dataset, dict_size, data, outlier, bcode, dims, ebs_L4, pred_err, comp_err);
	
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
     
    cout << "Best Performance: block_size = " << *blksz << ", vector_length = " << vecsz << ", sample_pct = " << sample_pct*100 << ", num_iters = " << num_iter << endl;  
	
	if (*blksz > 64 || *blksz < 8) *blksz = 64;
	return vecsz;
}


#endif
