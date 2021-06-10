// 200211

#ifndef WORKFLOW_HH
#define WORKFLOW_HH


#include "argument_parser/argparse.hh"
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

template <typename T, typename Q>
void Compress(argparse* ap,
              size_t&   num_outlier,
              bool      show_histo   = false) 
{
    std::string const& dataset   = ap->demo_dataset;
    std::string& finame          = ap->files.input_file;
    auto dims_L16                = InitializeDims(ap);
    auto eb_config               = new config_t(ap->dict_size, ap->eb);
    float sample_pct             = ap->sample_percentage;
    size_t len                   = dims_L16[LEN];
    int blksz                    = ap->block_size;
    int vecsz                    = ap->vector_length;
    int num_iterations           = ap->num_iterations;

    alignas(32) auto xdata   = new T[len]();
    alignas(32) auto outlier = new T[len]();
    alignas(32) auto code    = new Q[len]();

    if (ap->mode == "r2r") 
    {
        double value_range = GetDatumValueRange<float>(finame, dims_L16[LEN]);
        eb_config->ChangeToRelativeMode(value_range);
    }
    double const* const ebs_L4 = InitializeErrorBoundFamily(eb_config);
    auto tstart = hires::now(); // begin timing
    LogAll(log_info, "load", finame, len * (ap->dtype == "f32" ? sizeof(float) : sizeof(double)), "bytes,", ap->dtype);
    auto ldstart = hires::now(); 
    alignas(32) auto data     = io::ReadBinaryToNewArray<T>(finame, len);
    auto ldend   = hires::now(); 
    LogAll(log_dbg, "time loading datum:", static_cast<duration_t>(ldend - ldstart).count(), "sec");
    

    ////////////////////////////////////////////////////////////////////////////////
    // start of compression
    ////////////////////////////////////////////////////////////////////////////////
    double timing;
    if (ap->szwf.autotune)
    {
        LogAll(log_info, "begin autotune");
        auto astart = hires::now(); // begin timing
        vecsz = autotune_vector_len<T,Q>(ap, &blksz, &timing, data, outlier, code, dims_L16, ebs_L4);
        auto dims_L16 = InitializeDims(ap);
        auto aend = hires::now(); // begin timing
        LogAll(log_dbg, "autotune time:", static_cast<duration_t>(aend - astart).count(), "sec");
    }

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

    LogAll(log_info, "begin lossy-construction");
    auto pqstart = hires::now(); // begin timing
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

    auto pqend = hires::now(); //end timing
    double pqtime = static_cast<duration_t>(pqend - pqstart).count();
    double n_iterations = dims_L16[DIM0] * dims_L16[DIM1] * dims_L16[DIM2]; //perform op for each element in data
    double CGflops_s = ((n_iterations * CN_OPS) / pqtime) / pow(2,30); //n_iterations * n_ops / time
    double LGflops_s = ((n_iterations * LN_OPS) / pqtime) / pow(2,30); //n_iterations * n_ops / time
    LogAll(log_dbg, "pred+quant time:", pqtime, "sec");
    if (ap->verbose)
    {
        LogAll(log_dbg, "conservative gflops:", CGflops_s);
        LogAll(log_dbg, "leniant gflops:", LGflops_s);
    }

    // huffman encode
    uint8_t* out = NULL;
    size_t outSize = 0;
    LogAll(log_info, "begin huffman tree generation");

    auto hstart = hires::now(); // begin timing
    DV::HuffmanTree* tree = DV::createDefaultHuffmanTree();
    DV::encode_withTree(tree,code,len,&out,&outSize);
    auto hend = hires::now();
    double htime = static_cast<duration_t>(hend - hstart).count();
    if (ap->verbose) LogAll(log_dbg, "huffman time:", htime, "sec");

    auto tend = hires::now();
    LogAll(log_dbg, "compression time:", static_cast<duration_t>(tend - tstart).count(), "sec");

} // end Compression

template <typename T, typename Q>
void Decompress(std::string&        finame,  
                std::string const&  dataset,
                size_t const* const dims_L16,
                double const* const ebs_L4,
                size_t&             num_outlier,
                size_t              B,
                Q                   code,
                T                   xdata,
                T                   outlier,
                DV::HuffmanTree*    tree,
                uint8_t*            out,
	            int                 num_iterations,
	            float               sample_pct,
                bool                show_histo   = false) 
{            
    size_t len = dims_L16[LEN];
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
        Analysis::histogram(std::string("reconstructed datum"), xdata, len, 16);
    }

    io::WriteArrayToBinary(new string(finame + ".vecsz"), xdata, len);

} // end Decompress

}  // namespace interface

}  // namespace vecsz

#endif
