#pragma GCC option("arch=native","tune=native","no-zero-upper") //Enable AVX

#include <cstring>
#include <string>
#include <vector>

#include "argument_parser/argparse.hh"
#include "workflow.hh"
#include "types.hh"
#include "verify.hh"
#include "query.hh"
#include "dimensions.hh"

int main(int argc, char** argv) {
    std::string eb_mode, dataset, datum_path;
    double      eb;
    size_t      n_dim, d0, d1, d2, d3 = 1;
    int         num_iterations;
    float       sample_pct;
    size_t      nnz_outlier = 0;  // for calculating compression ratio

    // parse command line arguments
    auto ap = new ArgParse();
    ap->ParseVecszArgs(argc, argv);

    size_t blk = ap->block_size;
    const size_t DICT_SIZE = ap->dict_size;

    auto eb_config  = new config_t(DICT_SIZE, eb);
    auto dims_L16   = InitializeDims(ap);
    datum_path      = ap->files.input_file;
    dataset         = ap->demo_dataset;
    sample_pct      = ap->sample_percentage;
    num_iterations  = ap->num_iterations;
    eb_mode         = ap->mode;

    if (eb_mode == "r2r") {  // as of C++ 14, string is directly comparable?
        double value_range = GetDatumValueRange<float>(datum_path, dims_L16[LEN]);
        eb_config->ChangeToRelativeMode(value_range);
    }

    auto ebs_L4 = InitializeErrorBoundFamily(eb_config);


    if (ap->verbose) ap->PrintArgs();

#ifdef AUTOTUNE
    vecsz::interface::Compress<float, int>(datum_path, dataset, dims_L16, ebs_L4, nnz_outlier, blk, num_iterations, sample_pct, true);
#else
    vecsz::interface::Compress<float, int>(datum_path, dataset, dims_L16, ebs_L4, nnz_outlier, blk, true);
#endif
}
