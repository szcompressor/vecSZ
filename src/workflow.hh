// 200211

#ifndef WORKFLOW_HH
#define WORKFLOW_HH


#include "utils/io.hh"
#include "analysis.hh"
#include "dualquant.hh"
#include "huffman_cpu.hh"
#include "verify.hh"

#include "autotune.hh"



namespace pq  = vecsz::predictor_quantizer;
namespace DV = DesignVerification;

/* Used for timing */
using namespace std;
using namespace std::chrono;

namespace vecsz {

namespace interface {

namespace __loop {}

template <typename T, typename Q>
void Compress(std::string&        finame,  
              std::string const&  dataset,
              size_t const* const dims,
              double const* const ebs_L4,
              size_t&             num_outlier,
              size_t              B,
#ifdef AUTOTUNE
	          int                 num_iterations,
	          float               sample_pct,
#endif
              bool                show_histo   = false) 
{

    auto tstart = hires::now(); // begin timing
    size_t len = dims[LEN];

    alignas(32) auto data     = io::ReadBinaryToNewArray<T>(finame, len);
    alignas(32) auto data_cmp = io::ReadBinaryToNewArray<T>(finame, len);
    
    alignas(32) auto xdata   = new T[len]();
    alignas(32) auto outlier = new T[len]();
    alignas(32) auto code    = new Q[len]();

    ////////////////////////////////////////////////////////////////////////////////
    // start of compression
    ////////////////////////////////////////////////////////////////////////////////
    double timing;
    int vecsz = 256, blksz = B;
#ifdef AVX512
    vecsz = 512;
#endif
#ifdef AUTOTUNE
   astart = hires::now(); // begin timing
   vecsz = autotune_vector_len<T,Q,B>(num_iterations, sample_pct, &blksz, &timing, fine_massive, blocked, dataset, dims[CAP], data, outlier, code, dims, ebs_L4);
   auto dims_L16 = InitializeDemoDims(dataset, dims[CAP], blksz);
   aend = hires::now(); // begin timing
   double autotune_t = static_cast<duration_t>(aend - astart).count();
   cout << setprecision(6) << "Autotune Time: " << autotune_t << endl;
#else
    auto dims_L16 = dims;
#endif

    int CN_OPS, LN_OPS;
    if (dims_L16[nDIM] == 3) {
        CN_OPS = 15;
        LN_OPS = 18;
    }
    else if (dims_L16[nDIM] == 2) {
        CN_OPS = 11;
        LN_OPS = 14;
    }
    else {
        CN_OPS = 9;
        LN_OPS = 12;
    }

    auto cstart = hires::now(); // begin timing
    if (dims_L16[nDIM] == 1) 
    {
        #pragma omp parallel for
        for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
        {
            pq::c_lorenzo_1d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, b0, blksz, vecsz);
        }
    } 
    else if (dims_L16[nDIM] == 2) 
    {
        #pragma omp parallel for
        for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) 
        {
            for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
            {
                pq::c_lorenzo_2d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, b0, b1, blksz, vecsz);
            }
        }
    } 
    else if (dims_L16[nDIM] == 3) 
    {
        #pragma omp parallel for
        for (size_t b2 = 0; b2 < dims_L16[nBLK2]; b2++) 
        {
            for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) 
            {
                for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
                {
                    pq::c_lorenzo_3d1l<T, Q>(data, outlier, code, dims_L16, ebs_L4, b0, b1, b2, blksz, vecsz);
                }
            }
        }
    }

    auto cend = hires::now(); //end timing
    double tot_ctime = static_cast<duration_t>(cend - cstart).count();
    double n_iterations = dims_L16[DIM0] * dims_L16[DIM1] * dims_L16[DIM2]; //perform op for each element in data
    double CGflops_s = ((n_iterations * CN_OPS)/ tot_ctime) / pow(2,30); //n_iterations * n_ops / time
    double LGflops_s = ((n_iterations * LN_OPS)/ tot_ctime) / pow(2,30); //n_iterations * n_ops / time
//    retval = PAPI_hl_region_end("computation");
//    if ( retval != PAPI_OK ) {
//	cerr << "PAPI Error in region begin\n" << endl;
//	exit(1);
//   }
    cout << "Total PRED/QUANT time: " << tot_ctime << "s" << endl;
    cout << setprecision(6) << "Conservative GFLOPS: " << CGflops_s << endl;
    cout << setprecision(6) << "Leniant GFLOPS: " << LGflops_s << endl;

    //    io::write_binary_file(code, len, new string("/Users/jtian/WorkSpace/cuSZ/src/CLDMED.bincode"));

  // huffman encode
  uint8_t* out = NULL;
  size_t outSize = 0;
  cout << "Encoding prediction data..." << endl;

  auto hstart = hires::now(); // begin timing
  DV::HuffmanTree* tree = DV::createDefaultHuffmanTree();
  DV::encode_withTree(tree,code,len,&out,&outSize);
  auto hend = hires::now();
  double tot_huffman_time = static_cast<duration_t>(hend - hstart).count(); //end timing
  cout << "Total HUFFMAN time: " << tot_huffman_time << "s" << endl;

  //int huff_depth;
  //DV::hMD(tree,code,len,&huff_depth);

   // if (show_histo) {
   //     Analysis::histogram<int>(std::string("bincode/quant.code"), code, len, 8);
  //  }
 //   Analysis::getEntropy(code, len, 1024);

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
    auto tend = hires::now();
    double tot_sim_time = static_cast<duration_t>(tend - tstart).count(); //end timing
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

    if (dims_L16[nDIM] == 1) 
    {
        #pragma omp parallel for
        for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
        {
            pq::x_lorenzo_1d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0);
        }
    }
    else if (dims_L16[nDIM] == 2) 
    {
        #pragma omp parallel for
        for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) 
        {
            for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
            {
                pq::x_lorenzo_2d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1);
            }
        }
    } 
    else if (dims_L16[nDIM] == 3) 
    {
        #pragma omp parallel for
        for (size_t b2 = 0; b2 < dims_L16[nBLK2]; b2++) 
        {
            for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) 
            {
                for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
                {
                    pq::x_lorenzo_3d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1, b2);
                }
            }
        }
    }

    if (show_histo) {
        Analysis::histogram(std::string("original datum"), data_cmp, len, 16);
        Analysis::histogram(std::string("reconstructed datum"), xdata, len, 16);
    }

    cout << "\e[46mnum.outlier:\t" << num_outlier << "\e[0m" << endl;
    cout << setprecision(5) << "error bound: " << ebs_L4[EB] << endl;
	cout << setprecision(5) << "Compression Ratio: " << compressionRatio << endl;

    io::WriteBinaryFile(xdata, len, new string(finame + ".psz.cusz.out"));
    Analysis::VerifyData(xdata, data_cmp, len, 1);
*/
}

}  // namespace interface
}  // namespace vecsz

#endif
