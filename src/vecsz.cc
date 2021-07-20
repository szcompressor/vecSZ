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

int main(int argc, char** argv)
{
    std::string eb_mode, dataset, datum_path;
    double      eb;
    size_t      n_dim, d0, d1, d2, d3 = 1;
    int         num_iterations;
    float       sample_pct;
    size_t      nnz_outlier = 0;  // for calculating compression ratio

    // parse command line arguments
    auto ap = new ArgParse();
    ap->ParseVecszArgs(argc, argv);

    if (ap->verbose) GetMachineProperties();

    if (ap->szwf.lossy_compress) vecsz::interface::Compress<float, int>(ap, nnz_outlier);
    else if (ap->szwf.lossy_decompress) vecsz::interface::Decompress<float, int>(ap);
    else if (ap->szwf.lossy_dryrun) vecsz::interface::DryRun<float, int>(ap, nnz_outlier);
}
