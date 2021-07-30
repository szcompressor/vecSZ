#pragma GCC option("arch=native","tune=native","no-zero-upper") //Enable AVX

#include <cstring>
#include <string>
#include <vector>

#include "argument_parser/argparse.hh"
#include "vecsz_interface.hh"
#include "types.hh"
#include "dimensions.hh"

#include "utils/query.hh"
#include "utils/verify.hh"
#include "utils/io.hh"

int main(int argc, char** argv)
{
    // parse command line arguments
    auto ap = new ArgParse();
    ap->ParseVecszArgs(argc, argv);

    if (ap->verbose) GetMachineProperties();

    if (ap->szwf.lossy_compress)
    {
        size_t data_size = 0;
        size_t nnz_outlier = 0;
        unsigned long lossless_size = 0;

        InitializeDims(ap);
        // load data
        LogAll(log_info, "load", ap->files.input_file, ap->len * (ap->dtype == "f32" ? sizeof(float) : sizeof(double)), "bytes,", ap->dtype);
        auto ldstart = hires::now();
        alignas(32) auto data = io::ReadBinaryToNewArray<float>(ap->files.input_file, ap->len);
        auto ldend   = hires::now();
        LogAll(log_dbg, "time loading datum:", static_cast<duration_t>(ldend - ldstart).count(), "sec");

        // perform compression
        auto data_out = (unsigned char*)vecsz::interface::Compress<float, int>(ap, data, &nnz_outlier, &data_size, &lossless_size);

        // write data
        if (not ap->szwf.skip_write_output)
        {
            LogAll(log_info, "write", ap->files.output_file, data_size * sizeof(unsigned char), "bytes,", ap->dtype);
            io::WriteArrayToBinary<unsigned long>(ap->files.output_file, &(lossless_size), 1);
            io::AppendArrayToBinary<size_t>(ap->files.output_file, &(data_size), 1);
            io::AppendArrayToBinary<unsigned char>(ap->files.output_file, data_out, data_size);
        }
    }
    else if (ap->szwf.lossy_decompress)
    {
        // reading input file
        LogAll(log_info, "load", ap->files.input_file);
        auto lossless_size = io::ReadBinaryToNewArray<unsigned long>(ap->files.input_file, 1);
        auto target_size   = io::ReadBinaryToNewArrayPos<size_t>(ap->files.input_file, 1, sizeof(unsigned long));
        unsigned char* data_in;
        if (lossless_size[0] > 0) data_in = io::ReadBinaryToNewArrayPos<unsigned char>(ap->files.input_file, lossless_size[0], sizeof(size_t) + sizeof(unsigned long));
        else data_in = io::ReadBinaryToNewArrayPos<unsigned char>(ap->files.input_file, target_size[0], sizeof(size_t) + sizeof(unsigned long));

        //perform decompression
        auto data_out = (float*)vecsz::interface::Decompress<float,int>(ap, data_in, target_size[0], lossless_size[0]);

        if (not ap->szwf.skip_write_output)
        {
            LogAll(log_info, "write output to file", ap->files.output_file, ap->len * sizeof(float), "bytes");
            auto wstart = hires::now();
            io::WriteArrayToBinary<float>(ap->files.output_file, data_out, ap->len);
            auto wend = hires::now();
            if (ap->verbose) LogAll(log_dbg, "time writing datum:", static_cast<duration_t>(wend - wstart).count(), "sec");
        }

    }
    else if (ap->szwf.lossy_dryrun)
    {
        size_t data_size = 0;
        size_t nnz_outlier = 0;
        unsigned long lossless_size = 0;

        /************************************* BEGIN COMPRESSION ***********************************************/
        InitializeDims(ap);
        // load data
        LogAll(log_info, "load", ap->files.input_file, ap->len * (ap->dtype == "f32" ? sizeof(float) : sizeof(double)), "bytes,", ap->dtype);
        auto ldstart = hires::now();
        alignas(32) auto data_in  = io::ReadBinaryToNewArray<float>(ap->files.input_file, ap->len);
        alignas(32) auto ori_data = io::ReadBinaryToNewArray<float>(ap->files.input_file, ap->len);
        auto ldend   = hires::now();
        LogAll(log_dbg, "time loading datum:", static_cast<duration_t>(ldend - ldstart).count(), "sec");

        // perform compression
        auto cmp_data = (unsigned char*)vecsz::interface::Compress<float, int>(ap, data_in, &nnz_outlier, &data_size, &lossless_size);

        /*************************************** END COMPRESSION ***********************************************/

        /************************************* BEGIN DECOMPRESSION *********************************************/
        auto data_out = (float*)vecsz::interface::Decompress<float,int>(ap, cmp_data, data_size, lossless_size);

        if (not ap->szwf.skip_write_output)
        {
            LogAll(log_info, "write output to file", ap->files.output_file, ap->len * sizeof(float), "bytes");
            auto wstart = hires::now();
            io::WriteArrayToBinary<float>(ap->files.output_file, data_out, ap->len);
            auto wend = hires::now();
            if (ap->verbose) LogAll(log_dbg, "time writing datum:", static_cast<duration_t>(wend - wstart).count(), "sec");
        }

        delete[] cmp_data;
        /*************************************** END DECOMPRESSION *********************************************/

        /************************************* BEGIN ANALYSIS **************************************************/
        if (not ap->szwf.skip_verify)
        {
            if (ap->szwf.show_histo)
            {
                analysis::histogram(std::string("original datum"), ori_data, ap->len, 16);
                analysis::histogram(std::string("reconstructed data"), data_out, ap->len, 16);
            }
            analysis::VerifyData<float>(&(ap->stat), data_out, ori_data, ap->len);
            analysis::PrintMetrics<float>(&(ap->stat));
        }

        /*************************************** END ANALYSIS **************************************************/

        //clean up
        delete[] data_out;
        delete[] ori_data;
    }

    //if (ap->szwf.lossy_compress) vecsz::interface::Compress<float, int>(ap, nnz_outlier);
    //else if (ap->szwf.lossy_decompress) vecsz::interface::Decompress<float, int>(ap);
    //else if (ap->szwf.lossy_dryrun) vecsz::interface::DryRun<float, int>(ap, nnz_outlier);
}
