// 200211

#ifndef PSZ_WORKFLOW_HH
#define PSZ_WORKFLOW_HH

#include "io.hh"
#include "analysis.hh"
#include "psz_14.hh"
#include "psz_14blocked.hh"
#include "psz_dualquant_opt.hh"
#include "huffman_cpu.hh"
#include "verify.hh"
//#include <papi.h>

#include "autotune.hh"

#ifdef _3D
#define CN_OPS 15
#define LN_OPS 18
#elif _2D
#define CN_OPS 11
#define LN_OPS 14
#elif _1D
#define CN_OPS 9
#define LN_OPS 12
#endif


namespace PdQ  = pSZ::PredictionDualQuantization;
namespace PQRs = pSZ::PredictionQuantizationReconstructionSingleton;
namespace PQRb = pSZ::PredictionQuantizationReconstructionBlocked;
namespace DV = DesignVerification;

/* Used for timing */
using namespace std;
using namespace std::chrono;

high_resolution_clock::time_point now = high_resolution_clock::now();
high_resolution_clock::time_point timer = high_resolution_clock::now();
high_resolution_clock::time_point timer2 = high_resolution_clock::now();
#define TIME duration_cast<duration<double>>(high_resolution_clock::now() - now).count()
#define TIME2 duration_cast<duration<double>>(high_resolution_clock::now() - timer).count()
#define TIME3 duration_cast<duration<double>>(high_resolution_clock::now() - timer2).count()

namespace pSZ {

namespace FineMassiveSimulation {

namespace __loop {}

template <typename T, typename Q, size_t B>
void cx_sim(std::string&        finame,  //
            std::string const&  dataset,
            size_t const* const dims,
            double const* const ebs_L4,
            size_t&             num_outlier,
#ifdef AUTOTUNE
	    int			num_iterations,
	    float		sample_pct,
#endif
            bool                fine_massive = false,
            bool                blocked      = false,
            bool                show_histo   = false) {

    timer = high_resolution_clock::now(); // begin timing
    size_t len = dims[LEN];

    auto __attribute__((aligned(64))) data     = io::ReadBinaryFile<T>(finame, len);
    auto __attribute__((aligned(64))) data_cmp = io::ReadBinaryFile<T>(finame, len);

    T* pred_err = nullptr;
    T* comp_err = nullptr;
#ifdef PRED_COMP_ERR
    pred_err = new T[len]();
    comp_err = new T[len]();
#endif

    auto __attribute__((aligned(64))) xdata   = new T[len]();
    auto __attribute__((aligned(64))) outlier = new T[len]();
    auto __attribute__((aligned(64))) code    = new Q[len]();

    if (fine_massive)
        cout << "\e[46musing (blocked) dualquant\e[0m" << endl;
    else {
        cout << (blocked == true ? "\e[46musing blocked sz14\e[0m" : "\e[46musing non-blocked sz14\e[0m") << endl;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // start of compression
    ////////////////////////////////////////////////////////////////////////////////
    double timing;
    int vecsz = 256, blksz = B;
#ifdef AVX512
    vecsz = 512;
#endif
#ifdef AUTOTUNE
   timer2 = high_resolution_clock::now(); // begin timing
   vecsz = autotune_vector_len<T,Q,B>(num_iterations, sample_pct, &blksz, &timing, fine_massive, blocked, dataset, dims[CAP], data, outlier, code, dims, ebs_L4, pred_err, comp_err);
   auto dims_L16 = InitializeDemoDims(dataset, dims[CAP], blksz);
   double autotune_t = TIME3; 
   cout << setprecision(6) << "Autotune Time: " << autotune_t << endl;
#else
   auto dims_L16 = dims;
#endif

    // TODO omp version
    now = high_resolution_clock::now(); // begin timing
    const int ITERS = 1;
//    int retval = PAPI_hl_region_begin("computation");
//    if ( retval != PAPI_OK ) {
//	cerr << "PAPI Error in region begin\n" << endl;
//	exit(1);
//    }
for (int i = 0; i < ITERS; i++) { //FOR TESTING ONLY, REMOVE... make compression longer for profiling
    if (dims_L16[nDIM] == 1) {
        if (blocked) {
#pragma omp parallel for
            for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) {
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
            for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) {
                for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) {
                    if (fine_massive)
#ifdef REPBLK
                        PdQ::c_lorenzo_2d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, 0, 0, blksz, vecsz);
#else
                        PdQ::c_lorenzo_2d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0, b1, blksz, vecsz);
#endif

                    else
                        PQRb::c_lorenzo_2d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0, b1);
                }
            }

        } else {
            PQRs::c_lorenzo_2d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err);
        }
    } else if (dims_L16[nDIM] == 3) {
        if (blocked) {
#pragma omp parallel for
            for (size_t b2 = 0; b2 < dims_L16[nBLK2]; b2++) {
                for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) {
                    for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) {
                        if (fine_massive)
#ifdef REPBLK
                            PdQ::c_lorenzo_3d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, 0, 0, 0, blksz, vecsz);
#else
                            PdQ::c_lorenzo_3d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0, b1, b2, blksz, vecsz);
#endif
                        else
                            PQRb::c_lorenzo_3d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err, b0, b1, b2);
                    }
                }
            }
        } else {
            PQRs::c_lorenzo_3d1l<T, Q, B>(data, outlier, code, dims_L16, ebs_L4, pred_err, comp_err);
        }
    }
}
    double tot_cx_time = TIME; //end timing
    double n_iterations = dims_L16[DIM0] * dims_L16[DIM1] * dims_L16[DIM2]; //perform op for each element in data
    double CGflops_s = ((n_iterations * CN_OPS)/ tot_cx_time) / pow(2,30); //n_iterations * n_ops / time
    double LGflops_s = ((n_iterations * LN_OPS)/ tot_cx_time) / pow(2,30); //n_iterations * n_ops / time
//    retval = PAPI_hl_region_end("computation");
//    if ( retval != PAPI_OK ) {
//	cerr << "PAPI Error in region begin\n" << endl;
//	exit(1);
//   }
    cout << "Total PRED/QUANT time: " << tot_cx_time << "s" << endl;
    cout << "VEC AVG PRED/QUANT time: " << tot_cx_time/ITERS << "s" << endl;
    cout << setprecision(6) << "Conservative GFLOPS: " << CGflops_s/ITERS << endl;
    cout << setprecision(6) << "Leniant GFLOPS: " << LGflops_s/ITERS << endl;

    //    io::write_binary_file(code, len, new string("/Users/jtian/WorkSpace/cuSZ/src/CLDMED.bincode"));

  // huffman encode
  uint8_t* out = NULL;
  size_t outSize = 0;
  cout << "Encoding prediction data..." << endl;

  now = high_resolution_clock::now(); // begin timing
  DV::HuffmanTree* tree = DV::createDefaultHuffmanTree();
  DV::encode_withTree(tree,code,len,&out,&outSize);
  double tot_huffman_time = TIME; //end timing
  cout << "Total HUFFMAN time: " << tot_huffman_time << "s" << endl;
  cout << "Avg HUFFMAN time: " << tot_huffman_time/ITERS << "s" << endl;

  //int huff_depth;
  //DV::hMD(tree,code,len,&huff_depth);

   // if (show_histo) {
   //     Analysis::histogram<int>(std::string("bincode/quant.code"), code, len, 8);
  //  }
 //   Analysis::getEntropy(code, len, 1024);
#ifdef PRED_COMP_ERR
    Analysis::histogram<T>(std::string("pred.error"), pred_err, len, 8);  // TODO when changing to 8, seg fault
    Analysis::histogram<T>(std::string("comp.error"), comp_err, len, 16);
#endif
/*	FILE *fp2 = fopen("voutliers.csv","w");
	fprintf(fp2,"voutliers\n");
    for(size_t i = 0;i < len;i++) {
		if (outlier[i]) num_outlier += 1;
		fprintf(fp2, "%f\n",outlier[i]);
	}
	fclose(fp2);

	float pct_outlier = ((float)num_outlier/len) * 100;
	float compressionRatio = (len*sizeof(float))/(float)(outSize + num_outlier);
*/
    double tot_sim_time = TIME2; //end timing
	cout << "Total Runtime: " << tot_sim_time << "s" << endl;
//	cout << setprecision(5) << "Compression Ratio: " << compressionRatio << endl;
//	cout << "Number Outliers: " << num_outlier << endl;
//	cout << setprecision(10) << "Percentage Outliers: " << pct_outlier << endl;
//	cout << "Huffman Tree Depth: " << huff_depth << endl;

   //write to file
    //io::WriteBinaryFile(out, outSize, new string(finame + ".psz.tree.out"));
    //io::WriteBinaryFile(outlier, num_outlier * sizeof(T), new string(finame + ".psz.outlier.out"));

    ////////////////////////////////////////////////////////////////////////////////
    // start of decompression
    ////////////////////////////////////////////////////////////////////////////////
/*
	// huffman decode
	DesignVerification::decode_withTree(tree,out,len,code);

    if (dims_L16[nDIM] == 1) {
        if (blocked) {
#pragma omp parallel for
            for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) {
                if (fine_massive)
                    PdQ::x_lorenzo_1d1l<T, Q, B>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0);
                else
                    PQRb::x_lorenzo_1d1l<T, Q, B>(xdata, outlier, code, dims_L16, ebs_L4, b0);  // TODO __2EB
            }
        } else {
            PQRs::x_lorenzo_1d1l<T, Q, B>(xdata, outlier, code, dims_L16, ebs_L4);  // TODO __2EB
        }
    }
    if (dims_L16[nDIM] == 2) {
        if (blocked) {
#pragma omp parallel for
            for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) {
                for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) {
                    if (fine_massive)
                        PdQ::x_lorenzo_2d1l<T, Q, B>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1);
                    else
                        PQRb::x_lorenzo_2d1l<T, Q, B>(xdata, outlier, code, dims_L16, ebs_L4, b0, b1);  // TODO __2EB
                }
            }

        } else {
            PQRs::x_lorenzo_2d1l<T, Q, B>(xdata, outlier, code, dims_L16, ebs_L4);  // TODO __2EB
        }
    } else if (dims_L16[nDIM] == 3) {
        if (blocked) {
#pragma omp parallel for
            for (size_t b2 = 0; b2 < dims_L16[nBLK2]; b2++) {
                for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) {
                    for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) {
                        if (fine_massive)
                            PdQ::x_lorenzo_3d1l<T, Q, B>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1, b2);
                        else
                            PQRb::x_lorenzo_3d1l<T, Q, B>(xdata, outlier, code, dims_L16, ebs_L4, b0, b1, b2);  // TODO __2EB
                    }
                }
            }
        } else {
            PQRs::x_lorenzo_3d1l<T, Q, B>(xdata, outlier, code, dims_L16, ebs_L4);  // TODO __2EB
        }
    }


    if (show_histo) {
        Analysis::histogram(std::string("original datum"), data_cmp, len, 16);
        Analysis::histogram(std::string("reconstructed datum"), xdata, len, 16);
    }

    cout << "\e[46mnum.outlier:\t" << num_outlier << "\e[0m" << endl;
    cout << setprecision(5) << "error bound: " << ebs_L4[EB] << endl;
	cout << setprecision(5) << "Compression Ratio: " << compressionRatio << endl;

    if (fine_massive) {
        io::WriteBinaryFile(xdata, len, new string(finame + ".psz.cusz.out"));
    } else if (blocked == true and fine_massive == false) {
        io::WriteBinaryFile(xdata, len, new string(finame + ".psz.sz14blocked.out"));
    } else if (blocked == false and fine_massive == false) {
        io::WriteBinaryFile(xdata, len, new string(finame + ".psz.sz14.out"));
        io::WriteBinaryFile(pred_err, len, new string(finame + ".psz.sz14.prederr"));
        io::WriteBinaryFile(comp_err, len, new string(finame + ".psz.sz14.xerr"));
    }
    Analysis::VerifyData(xdata, data_cmp, len, 1);
*/
}

}  // namespace FineMassiveSimulation
}  // namespace pSZ

#endif
