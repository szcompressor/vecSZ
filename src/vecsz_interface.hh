// 200211

#ifndef WORKFLOW_HH
#define WORKFLOW_HH


#include "argument_parser/argparse.hh"
#include "autotune.hh"
#include "dualquant.hh"
#include "huffman.hh"
#include "lossless.hh"
#include "utils/analysis.hh"
#include "utils/verify.hh"
#include "utils/io.hh"
#include "utils/datapack.hh"
#include "utils/padding.hh"


namespace pq  = vecsz::predictor_quantizer;
namespace DV = DesignVerification;

/* Used for timing */
using namespace std;
using namespace std::chrono;

namespace vecsz {

namespace interface {

template <typename T, typename Q>
void* Compress(argparse* ap,
               T* data_in,
               size_t* num_outlier,
               size_t* data_size,
               unsigned long* lossless_size)
{
    auto tstart              = hires::now();
    std::string& finame      = ap->files.input_file;
    auto dims_L16            = InitializeDims(ap);
    auto eb_config           = new config_t(ap->dict_size, ap->eb);
    size_t len               = dims_L16[LEN];
    int blksz                = ap->block_size;
    int vecsz                = ap->vector_length;

    alignas(32) auto outlier = new T[len]();
    alignas(32) auto code    = new Q[len]();

    // adjust error bound to relative if necessary
    if (ap->mode == "r2r") { double value_range = GetDatumValueRange<float>(finame, dims_L16[LEN]);
        eb_config->ChangeToRelativeMode(value_range);
    }
    double const* const ebs_L4 = InitializeErrorBoundFamily(eb_config);

    if (ap->verbose)
    {
        ap->PrintArgs();
    }

    ////////////////////////////////////////////////////////////////////////////////
    // start of compression
    ////////////////////////////////////////////////////////////////////////////////
    double timing;
    if (ap->szwf.autotune)
    {
        LogAll(log_info, "begin autotune");
        auto astart = hires::now(); // begin timing
        vecsz = autotune_vector_len<T,Q>(ap, &blksz, &timing, data_in, outlier, code, dims_L16, ebs_L4);
        dims_L16 = InitializeDims(ap);
        auto aend = hires::now(); // begin timing
        LogAll(log_dbg, "autotune time:", static_cast<duration_t>(aend - astart).count(), "sec");
    }

    // find global padding value using entire dataset
    size_t pad_idx = 0;
    alignas(32) T* pad_vals = NULL;
    if (ap->szwf.global_padding)
    {
        size_t pad_dims[3] {dims_L16[DIM0], dims_L16[DIM1], dims_L16[DIM2]};
        ap->pad_constant = padding::find_pad_value<T>(data_in, ap->pad_type, pad_dims, ap->pad_constant);
        ap->pad_type = CONST_VAL;
        pad_vals = new T[1];
        pad_vals[pad_idx++] = round(ap->pad_constant * ebs_L4[EBx2_r]); // account for prequant
    }
    else if (ap->szwf.block_padding)
    {
        auto pad_vals_arr = new T[ap->nblk4._0 * ap->nblk4._1 * ap->nblk4._2 * ap->nblk4._3];
        pad_vals = pad_vals_arr;
    }
    else if (ap->szwf.edge_padding)
    {
        auto pad_vals_arr = new T[ap->ndim * ap->nblk4._0 * ap->nblk4._1 * ap->nblk4._2 * ap->nblk4._3];
        pad_vals = pad_vals_arr;
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
            pq::c_lorenzo_1d1l<T, Q>(data_in, outlier, code, dims_L16, ebs_L4, b0, blksz, vecsz);
        }
    }
    else if (dims_L16[nDIM] == 2)
    {
        #pragma omp parallel for
        for (size_t b1 = 0; b1 < dims_L16[nBLK1]; b1++)
        {
            for (size_t b0 = 0; b0 < dims_L16[nBLK0]; b0++)
            {
                pq::c_lorenzo_2d1l<T, Q>(data_in, outlier, code, dims_L16, ebs_L4, b0, b1, blksz, vecsz, ap->szwf, ap->pad_constant, ap->pad_type, pad_vals, &pad_idx);
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
                    pq::c_lorenzo_3d1l<T, Q>(data_in, outlier, code, dims_L16, ebs_L4, b0, b1, b2, blksz, vecsz, ap->szwf, ap->pad_constant, ap->pad_type, pad_vals, &pad_idx);
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
        if (code[i] == 0) outliers[(*num_outlier)++] = outlier[i];
    }
    if (ap->verbose) LogAll(log_dbg, "number of outliers:", *num_outlier);

    auto tend = hires::now();
    LogAll(log_info, "complete lossy compression:");
    LogAll(log_dbg, "compression time:", static_cast<duration_t>(tend - tstart).count(), "sec");

    //if (ap->verbose)
    //{
	    analysis::get_outliers(dims_L16, code, ap->block_size, *num_outlier);
    //}

    if (ap->szwf.show_histo)
    {
        analysis::histogram(std::string("bincode/quant.code"), code, len, 8);
        analysis::getEntropy<Q>(code,len,ap->dict_size);
    }

        // organize necessary metadata and write to file
	    *data_size = (pad_idx + *num_outlier) * sizeof(T) + huff_size * sizeof(uint8_t) + 1 * sizeof(double); /* OUTLIER + HUFF_ENC_DATA + EB */

	    //dimension info size
        auto dims_metadata = new int[dims_L16[nDIM] + 1];
        dims_metadata[0] = dims_L16[nDIM];
        for (size_t i = 0; i < dims_L16[nDIM]; i++) dims_metadata[i + 1] = dims_L16[i];
	    *data_size += (dims_L16[nDIM] + 1) * sizeof(int); /* + DIMS */

	    // //length of compressed data size
        // auto len_metadata = new size_t[2];
        // len_metadata[0] = huff_size;
        // len_metadata[1] = *num_outlier;
	    // *data_size += 2 * sizeof(size_t); /* + LEN_METADATA */

	    //length of compressed data size
        auto len_metadata = new size_t[3];
        len_metadata[0] = pad_idx;
        len_metadata[1] = huff_size;
        len_metadata[2] = *num_outlier;
	    *data_size += 3 * sizeof(size_t); /* + LEN_METADATA */

	    auto data_out = datapack::pack(*data_size, dims_metadata, dims_L16[nDIM], len_metadata, 3, ap->eb, pad_vals, huff_result, outliers);

	    // lossless pass
	    LogAll(log_dbg, "compression ratio:", (float)(len * sizeof(T)) / *data_size);
        if (ap->szwf.lossless_pass)
        {
	        unsigned char* lossless_out = NULL;
            *lossless_size = 0;

	        if (ap->szwf.lossless_gzip) *lossless_size = vecsz::lossless::sz_lossless_compress(GZIP_COMPRESSOR, 3, data_out, *data_size, &lossless_out);
            else *lossless_size = vecsz::lossless::sz_lossless_compress(ZSTD_COMPRESSOR, 3, data_out, *data_size, &lossless_out);
	        LogAll(log_dbg, "compression ratio after lossless pass:", (float)(len * sizeof(T)) / *lossless_size);

            return lossless_out;
        }

        //clean-up
        delete[] code;
        delete[] data_in;
        delete[] outlier;
        delete[] outliers;
        delete[] len_metadata;
        delete[] dims_metadata;

        return data_out;
}

template <typename T, typename Q>
void* Decompress(argparse*      ap,
                 unsigned char* tmp_data,
                 size_t         target_size,
                 unsigned long  lossless_size)
{
    // perform lossless pass
    unsigned char *data_in = NULL;
    int lossless_pass = vecsz::lossless::is_lossless_compressed_data(tmp_data, lossless_size);
    if (lossless_pass != -1)
    {
        LogAll(log_info, "initiate lossless pass");
        vecsz::lossless::sz_lossless_decompress(lossless_pass, tmp_data, lossless_size, &data_in, target_size);
    }
    else
    {
        data_in = tmp_data;
    }

    //assign values
    int      ndims;
    double   eb;
    int*     dims;
    size_t   pad_length;
    uint8_t* huffman;
    T*       outliers;
    T*       pad_vals;

    LogAll(log_info, "unpack input data");
    datapack::unpack<uint8_t,T>(data_in, &ndims, &dims, &eb, &pad_vals, &pad_length, &huffman, &outliers);

    ap->ndim = ndims;
    if (ap->ndim == 1) {ap->dim4._0 = dims[0]; ap->dim4._1 = 1; ap->dim4._2 = 1; ap->dim4._3 = 1;}
    if (ap->ndim == 2) {ap->dim4._0 = dims[0]; ap->dim4._1 = dims[1]; ap->dim4._2 = 1; ap->dim4._3 = 1;}
    if (ap->ndim == 3) {ap->dim4._0 = dims[0]; ap->dim4._1 = dims[1]; ap->dim4._2 = dims[2]; ap->dim4._3 = 1;}
    if (ap->ndim == 4) {ap->dim4._0 = dims[0]; ap->dim4._1 = dims[1]; ap->dim4._2 = dims[2]; ap->dim4._3 = dims[3];}

    ap->eb = eb;

    size_t const* const dims_L16       = InitializeDims(ap);
    auto                eb_config      = new config_t(ap->dict_size, ap->eb);
    size_t              blksz          = ap->block_size;
    std::string&        finame         = ap->files.input_file;
    size_t              len            = dims_L16[LEN];
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

    // determine padding configuration
    size_t pad_idx = 0;
    size_t num_blocks = ap->nblk4._0 * ap->nblk4._1 * ap->nblk4._2 * ap->nblk4._3;
    if (pad_length == 1) ap->szwf.global_padding = true;
    else if (pad_length == num_blocks) ap->szwf.block_padding = true;
    else if (pad_length == ap->ndim * num_blocks) ap->szwf.edge_padding = true;
    else
    {
        ap->szwf.global_padding = false;
        ap->szwf.block_padding  = false;
        ap->szwf.edge_padding   = false;
    }

    cout << log_info << "using " << std::string((ap->szwf.global_padding)?"global":
                                    ((ap->szwf.block_padding)?"block":
                                    ((ap->szwf.edge_padding)?"edge":"no"))) << " padding" << endl;

	// huffman decode
    LogAll(log_info, "invoke huffman tree reconstruction");
    auto hstart = hires::now();
    tree = DV::createHuffmanTree(ap->dict_size * 2);
	DV::decode_withTree(tree,huffman,len,code);
    auto hend = hires::now();
    if (ap->verbose) LogAll(log_dbg, "huffman reconstruction time:", static_cast<duration_t>(hend - hstart).count(), "sec");

    // reconstruct outliers array
    size_t oi = 0;
    for (size_t i = 0; i < len; i++)
    {
       outlier[i] = (code[i] == 0) ? outliers[oi++] : 0;
    }
    if (ap->verbose) LogAll(log_dbg, "number of outliers:", oi);

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
                pq::x_lorenzo_2d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1, blksz, ap->szwf, &pad_idx, pad_vals);
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
                    pq::x_lorenzo_3d1l<T, Q>(xdata, outlier, code, dims_L16, ebs_L4[EBx2], b0, b1, b2, blksz, ap->szwf, &pad_idx, pad_vals);
                }
            }
        }
    }
    auto pqend = hires::now();
    LogAll(log_dbg, "reverse pred+quant time:", static_cast<duration_t>(pqend - pqstart).count(), "sec");


   // display histogramming info for reconstructed data
   if (ap->szwf.show_histo)
   {
       analysis::histogram(std::string("bincode/quant.code"), code, len, 8);
       analysis::histogram(std::string("reconstructed datum"), xdata, len, 16);
    }

    // clean up
    delete[] code;
    delete[] dims;
    delete[] huffman;
    delete[] outlier;
    delete[] outliers;

    return xdata;
}

}  // namespace interface

}  // namespace vecsz

#endif
