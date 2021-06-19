// 200211

#ifndef WORKFLOW_HH
#define WORKFLOW_HH


#include "argument_parser/argparse.hh"
#include "autotune.hh"
#include "dualquant.hh"
#include "huffman.hh"
#include "utils/analysis.hh"
#include "utils/verify.hh"
#include "utils/io.hh"



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

    alignas(32) auto outlier = new T[len]();
    alignas(32) auto code    = new Q[len]();

    if (ap->mode == "r2r") 
    {
        double value_range = GetDatumValueRange<float>(finame, dims_L16[LEN]);
        eb_config->ChangeToRelativeMode(value_range);
    }
    double const* const ebs_L4 = InitializeErrorBoundFamily(eb_config);

    if (ap->verbose)
    {
        ap->PrintArgs();
    }
    
    auto tstart  = hires::now(); // begin timing

    // load data
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
    uint8_t* huff_result = NULL;
    size_t huff_size = 0;
    LogAll(log_info, "begin huffman tree generation");

    auto hstart = hires::now(); // begin timing
    DV::HuffmanTree* tree = DV::createHuffmanTree(ap->dict_size * 2);
    DV::encode_withTree(tree,code,len,&huff_result,&huff_size);
    auto hend = hires::now();
    double htime = static_cast<duration_t>(hend - hstart).count();
    if (ap->verbose) LogAll(log_dbg, "huffman time:", htime, "sec");

    for (size_t i = 0; i < len; i++)
    {
        data[num_outlier] = outlier[i];
        num_outlier += (outlier[i] == 0) ? 0 : 1;
    }
    if (ap->verbose) LogAll(log_dbg, "number of outliers:", num_outlier);

    auto tend = hires::now();
    LogAll(log_info, "complete lossy compression:");
    LogAll(log_dbg, "compression time:", static_cast<duration_t>(tend - tstart).count(), "sec");

    LogAll(log_info, "write", ap->files.output_file, len * (ap->dtype == "f32" ? sizeof(float) : sizeof(double)), "bytes,", ap->dtype);
    auto wstart = hires::now();

    { // organize necessary metadata and write to file
        auto dims_metadata = new int[dims_L16[nDIM] + 1];
        dims_metadata[0] = dims_L16[nDIM];
        for (int i = 0; i < dims_L16[nDIM]; i++) dims_metadata[i + 1] = dims_L16[i];

        auto len_metadata = new size_t[2];
        len_metadata[0] = huff_size;
        len_metadata[1] = num_outlier;

        io::WriteArrayToBinary<int>(ap->files.output_file, dims_metadata, dims_L16[nDIM] + 1); 
        io::AppendArrayToBinary<double>(ap->files.output_file, &(ap->eb), 1);
        io::AppendArrayToBinary<size_t>(ap->files.output_file, len_metadata, 2); 
        io::AppendArrayToBinary<uint8_t>(ap->files.output_file, huff_result, huff_size);
        io::AppendArrayToBinary<T>(ap->files.output_file, outlier, num_outlier);
    }

    auto wend = hires::now();
    LogAll(log_dbg, "time writing datum:", static_cast<duration_t>(wend - wstart).count(), "sec");
    
} // end Compression

template <typename T, typename Q>
void Decompress(argparse*        ap,
                bool             show_histo = false)
{            
    auto dstart = hires::now();

    uint8_t* huffman;

    auto ldstart = hires::now();
    LogAll(log_info, "load", ap->files.input_file);
    // reading input file
    auto ndims    = io::ReadBinaryToNewArray<int>(ap->files.input_file, 1); 
    auto dims     = io::ReadBinaryToNewArrayPos<int>(ap->files.input_file, ndims[0], sizeof(int)); 
    auto eb       = io::ReadBinaryToNewArrayPos<double>(ap->files.input_file, 1, (ndims[0] + 1) * sizeof(int));
    auto lengths  = io::ReadBinaryToNewArrayPos<size_t>(ap->files.input_file, 2, sizeof(double) + (ndims[0] + 1) * sizeof(int));
    huffman       = io::ReadBinaryToNewArrayPos<uint8_t>(ap->files.input_file, lengths[0], 2 * sizeof(size_t) + sizeof(double) + (ndims[0] + 1) * sizeof(int));
    auto outliers = io::ReadBinaryToNewArrayPos<float>(ap->files.input_file, lengths[1], lengths[0] * sizeof(uint8_t) + 2 * sizeof(size_t) + sizeof(double) + (ndims[0] + 1) * sizeof(int));
    auto ldend    = hires::now();
    LogAll(log_dbg, "time loading datum:", static_cast<duration_t>(ldend - ldstart).count(), "sec");

    ap->ndim = ndims[0];
    if (ap->ndim == 1) {ap->dim4._0 = dims[0]; ap->dim4._1 = 1; ap->dim4._2 = 1; ap->dim4._3 = 1;}
    if (ap->ndim == 2) {ap->dim4._0 = dims[0]; ap->dim4._1 = dims[1]; ap->dim4._2 = 1; ap->dim4._3 = 1;}
    if (ap->ndim == 3) {ap->dim4._0 = dims[0]; ap->dim4._1 = dims[1]; ap->dim4._2 = dims[2]; ap->dim4._3 = 1;}
    if (ap->ndim == 4) {ap->dim4._0 = dims[0]; ap->dim4._1 = dims[1]; ap->dim4._2 = dims[2]; ap->dim4._3 = dims[3];}

    ap->eb = eb[0];
    
    size_t const* const dims_L16       = InitializeDims(ap);
    auto                eb_config      = new config_t(ap->dict_size, ap->eb);
    size_t              blksz          = ap->block_size;
    std::string&        finame         = ap->files.input_file;  
    size_t              len            = dims_L16[LEN];
    int                 num_iterations = ap->num_iterations;
    float               sample_pct     = ap->sample_percentage;
    DV::HuffmanTree*    tree;
    auto                xdata          = new T[len];
    auto                outlier        = new T[len];
    auto                code           = new Q[len];

    if (ap->mode == "r2r") 
    {
        double value_range = GetDatumValueRange<float>(finame, dims_L16[LEN]);
        eb_config->ChangeToRelativeMode(value_range);
    }
    double const* const ebs_L4 = InitializeErrorBoundFamily(eb_config);

	// huffman decode
    LogAll(log_info, "invoke huffman tree reconstruction");
    auto hstart = hires::now();
    tree = DV::createHuffmanTree(ap->dict_size * 2);
	DV::decode_withTree(tree,huffman,len,code);
    auto hend = hires::now();
    if (ap->verbose) LogAll(log_dbg, "huffman reconstruction time:", static_cast<duration_t>(hend - hstart).count(), "sec");

    //reconstruct outliers array
    for (size_t i = 0, oi = 0; i < len; i++)
    {
       outlier[i] = (code[i] == 0) ? outliers[oi++] : 0;
    }


    LogAll(log_info, "invoke reverse prediction-quantization");
    auto pqstart = hires::now();
    // begin reverse prediction-quantization
    if (dims_L16[nDIM] == 1) 
    {
        #pragma omp parallel for
        for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
        {
            pq::x_lorenzo_1d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, blksz);
        }
    }
    else if (dims_L16[nDIM] == 2) 
    {
        #pragma omp parallel for
        for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) 
        {
            for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
            {
                pq::x_lorenzo_2d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1, blksz);
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
                    pq::x_lorenzo_3d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1, b2, blksz);
                }
            }
        }
    }
    auto pqend = hires::now();
    LogAll(log_dbg, "reverse pred+quant time:", static_cast<duration_t>(pqend - pqstart).count(), "sec");

    LogAll(log_info, "write output to file",ap->files.output_file,len,"bytes");
    auto wstart = hires::now();
    io::WriteArrayToBinary(ap->files.output_file, xdata, len);
    auto wend = hires::now();
    if (ap->verbose) LogAll(log_dbg, "time writing datum:", static_cast<duration_t>(wend - wstart).count(), "sec");

    auto dend = hires::now();
    LogAll(log_dbg, "decompression time:", static_cast<duration_t>(dend - dstart).count(), "sec");

} // end Decompress

template <typename T, typename Q>
void DryRun(argparse* ap,
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

    alignas(32) auto outlier = new T[len]();
    alignas(32) auto code    = new Q[len]();

    if (ap->mode == "r2r") 
    {
        double value_range = GetDatumValueRange<float>(finame, dims_L16[LEN]);
        eb_config->ChangeToRelativeMode(value_range);
    }
    double const* const ebs_L4 = InitializeErrorBoundFamily(eb_config);

    if (ap->verbose)
    {
        ap->PrintArgs();
    }
    
    auto tstart  = hires::now(); // begin timing

    // load data
    LogAll(log_info, "load", finame, len * (ap->dtype == "f32" ? sizeof(float) : sizeof(double)), "bytes,", ap->dtype);
    auto ldstart = hires::now(); 
    alignas(32) auto data     = io::ReadBinaryToNewArray<T>(finame, len);
    auto ldend   = hires::now(); 
    LogAll(log_dbg, "time loading datum:", static_cast<duration_t>(ldend - ldstart).count(), "sec");
    

    // -------------------------------------BEGIN-COMPRESSION----------------------------------------
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
    uint8_t* huff_result = NULL;
    size_t huff_size = 0;
    LogAll(log_info, "begin huffman tree generation");

    auto hstart = hires::now(); // begin timing
    DV::HuffmanTree* tree = DV::createHuffmanTree(ap->dict_size * 2);
    DV::encode_withTree(tree,code,len,&huff_result,&huff_size);
    auto hend = hires::now();
    double htime = static_cast<duration_t>(hend - hstart).count();
    if (ap->verbose) LogAll(log_dbg, "huffman time:", htime, "sec");

    auto outliers = new T[len];
    for (size_t i = 0; i < len; i++)
    {
        outliers[num_outlier] = outlier[i];
        num_outlier += (outlier[i] == 0) ? 0 : 1;
    }
    if (ap->verbose) LogAll(log_dbg, "number of outliers:", num_outlier);

    auto tend = hires::now();
    LogAll(log_info, "complete lossy compression:");
    LogAll(log_dbg, "compression time:", static_cast<duration_t>(tend - tstart).count(), "sec");

    delete[] outlier;
    delete[] code;

    // ---------------------------------------END-COMPRESSION----------------------------------------

    // -------------------------------------BEGIN-DECOMPRESSION--------------------------------------
    auto dstart  = hires::now();
    auto xdata   = new T[len];
    code    = new Q[len];
    outlier = new T[len];


	// huffman decode
    LogAll(log_info, "invoke huffman tree reconstruction");
    hstart = hires::now();
    tree = DV::createHuffmanTree(ap->dict_size * 2);
    DV::decode_withTree(tree,huff_result,len,code);
    hend = hires::now();
    if (ap->verbose) LogAll(log_dbg, "huffman reconstruction time:", static_cast<duration_t>(hend - hstart).count(), "sec");

    //reconstruct outliers array
    for (size_t i = 0, oi = 0; i < len; i++)
    {
       outlier[i] = (code[i] == 0) ? outliers[oi++] : 0;
    }

    LogAll(log_info, "invoke reverse prediction-quantization");
    pqstart = hires::now();
    // begin reverse prediction-quantization
    if (dims_L16[nDIM] == 1) 
    {
        #pragma omp parallel for
        for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
        {
            pq::x_lorenzo_1d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, blksz);
        }
    }
    else if (dims_L16[nDIM] == 2) 
    {
        #pragma omp parallel for
        for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++) 
        {
            for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++) 
            {
                pq::x_lorenzo_2d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1, blksz);
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
                    pq::x_lorenzo_3d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1, b2, blksz);
                }
            }
        }
    }
    pqend = hires::now();
    LogAll(log_dbg, "reverse pred+quant time:", static_cast<duration_t>(pqend - pqstart).count(), "sec");

    auto dend = hires::now();
    LogAll(log_dbg, "decompression time:", static_cast<duration_t>(dend - dstart).count(), "sec");

    // --------------------------------------END-DECOMPRESSION---------------------------------------
    
    // -------------------------------------BEGIN-VERIFICATION---------------------------------------
    if (show_histo)
    {
        analysis::histogram(std::string("original datum"), data, len, 16);
        analysis::histogram(std::string("reconstructed datum"), xdata, len, 16);
    }
    analysis::VerifyData<T>(&(ap->stat), xdata, data, len);
    analysis::PrintMetrics<T>(&(ap->stat));
    // ---------------------------------------END-VERIFICATION---------------------------------------
}

}  // namespace interface

}  // namespace vecsz

#endif
