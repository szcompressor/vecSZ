#pragma GCC option("arch=native","tune=native","no-zero-upper") //Enable AVX

#include <cstring>
#include <string>
#include <vector>

#include "SDRB.hh"
#include "psz_workflow_opt.hh"
#include "types.hh"
#include "verify.hh"


namespace fm = pSZ::FineMassiveSimulation;

// const size_t DICT_SIZE = 1024;
const size_t DICT_SIZE = 4096;

int main(int argc, char** argv) {
    std::string eb_mode, dataset, datum_path;
    bool        if_blocking, if_dualquant;
    double      mantissa, exponent;
    size_t      n_dim, d0, d1, d2, d3 = 1;
#ifdef AUTOTUNE
    int 	num_iterations;
    float	sample_pct;
#endif

#if defined(_1D)
    cout << "\e[46mThis program is working for 1D datasets.\e[0m" << endl;
#elif defined(_2D)
    cout << "\e[46mThis program is working for 2D datasets.\e[0m" << endl;
#elif defined(_3D)
    cout << "\e[46mThis program is working for 3D datasets.\e[0m" << endl;
#endif

#ifdef AUTOTUNE
    if (argc != 8 && argc != 11 && argc != 10) {
#else
    if (argc != 8 && argc != 11) {
#endif
        cout << "./<program> <abs|rel2range OR r2r> <mantissa> <exponent> <if blocking> <if dualquant> <dataset> <datum_path path>" << endl;
        cout << "OR" << endl;
        cout << "./<program> <abs|rel2range OR r2r> <mantissa> <exponent> <if blocking> <if dualquant> <-1/-2/-3> nx [ny [nz]] <datum_path path>" << endl;
        cout << "supported dimension and datasets" << endl;
        cout << "\t1D\t./psz1d r2r 1.23 -4.56 <noblk|yesblk> <nodq|dq> <hacc> /path/to/vx.f32" << endl;
        cout << "\t2D\t./psz2d r2r 1.23 -4.56 <noblk|yesblk> <nodq|dq> <cesm> /path/to/CLDHGH_1_1800_3600.f32" << endl;
        cout << "\t3D\t./psz3d r2r 1.23 -4.56 <noblk|yesblk> <nodq|dq> <hurricane|nyx|qmc|qmcpre> "
                "/path/to/CLOUDf48.bin.f32"
             << endl;
        exit(0);
    } else {
        eb_mode      = std::string(argv[1]);
        mantissa     = std::stod(argv[2]);
        exponent     = std::stod(argv[3]);
        if_blocking  = std::string(argv[4]) == "yesblk";
        if_dualquant = std::string(argv[5]) == "dq";
#ifdef AUTOTUNE
	if (argc == 8) {
            dataset        = std::string(argv[6]);
            datum_path     = std::string(argv[7]);
	    num_iterations = 10;
	    sample_pct	   = (float) 1 / 100;
	}
	if (argc == 10) {
            dataset        = std::string(argv[6]);
            datum_path     = std::string(argv[7]);
	    num_iterations = std::stoi(argv[8]);
	    sample_pct     = std::stof(argv[9]) / 100;
	}
#else
	if (argc == 8) {
            dataset = std::string(argv[6]);
            datum_path   = std::string(argv[7]);
        }
#endif
        else {
            n_dim        = atoi(&argv[6][1]);
            d0           = atoi(argv[7]);
            d1           = atoi(argv[8]);
            d2           = atoi(argv[9]);
            datum_path   = std::string(argv[10]);
        }
    }
    //for_each(argv, argv + 8, [](char* i) { cout << i << " "; });
    for_each(argv, argv + 8, [](auto i) { cout << i << " "; });
    cout << endl;
    size_t blk = 64;
    auto eb_config  = new config_t(DICT_SIZE, mantissa, exponent);
#ifdef AUTOTUNE
    auto dims_L16 = (argc == 8 || argc == 10) ? InitializeDemoDims(dataset, DICT_SIZE,blk):InitializeDims(DICT_SIZE,n_dim,d0,d1,d2,d3,blk);
#else 
    auto dims_L16 = (argc == 8) ? InitializeDemoDims(dataset, DICT_SIZE,blk):InitializeDims(DICT_SIZE,n_dim,d0,d1,d2,d3,blk);
#endif
    //auto dims_L16 = (argc == 8) ? InitializeDemoDims(dataset, DICT_SIZE):InitializeDims(DICT_SIZE,n_dim,d0,d1,d2,d3);
    printf("%-20s%s\n", "filename", datum_path.c_str());
    printf("%-20s%lu\n", "filesize", dims_L16[LEN] * sizeof(float));
    if (eb_mode == "r2r") {  // as of C++ 14, string is directly comparable?
        double value_range = GetDatumValueRange<float>(datum_path, dims_L16[LEN]);
        eb_config->ChangeToRelativeMode(value_range);
    }
    eb_config->debug();
    //    size_t c_byteSize;
    size_t num_outlier = 0;  // for calculating compression ratio

    cout << "block size:\t" << blk << endl;
    auto ebs_L4 = InitializeErrorBoundFamily(eb_config);

#ifdef AUTOTUNE
    fm::cx_sim<float, int>(datum_path, dataset, dims_L16, ebs_L4, num_outlier, blk, num_iterations, sample_pct, if_dualquant, if_blocking, true);
#else
    fm::cx_sim<float, int>(datum_path, dataset, dims_L16, ebs_L4, num_outlier, blk, if_dualquant, if_blocking, true);
#endif
}
