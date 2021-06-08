#include <cstring>
#include <string>
#include <vector>

#include "SDRB.hh"
#include "psz_workflow.hh"
#include "types.hh"
#include "verify.hh"

namespace fm = pSZ::FineMassiveSimulation;

// const size_t DICT_SIZE = 1024;
const size_t DICT_SIZE = 4096;
#if defined(_1D)
#ifdef B4
const int BLK = 4;
#elif B6
const int BLK = 6;
#elif B8
const int BLK = 8;
#elif B16
const int BLK = 16;
#elif B32
const int BLK = 32;
#elif B64
const int BLK = 64;
#elif B128
const int BLK = 128;
#elif B256
const int BLK = 256;
#endif
#elif defined(_2D)
#ifdef B4
const int BLK = 4;
#elif B6
const int BLK = 6;
#elif B8
const int BLK = 8;
#elif B16
const int BLK = 16;
#elif B32
const int BLK = 32;
#elif B64
const int BLK = 64;
#elif B128
const int BLK = 128;
#elif B256
const int BLK = 256;
#endif
#elif defined(_3D)
#ifdef B4
const int BLK = 4;
#elif B6
const int BLK = 6;
#elif B8
const int BLK = 8;
#elif B16
const int BLK = 16;
#elif B32
const int BLK = 32;
#elif B64
const int BLK = 64;
#elif B128
const int BLK = 128;
#elif B256
const int BLK = 256;
#endif
#endif

int main(int argc, char** argv) {
    std::string eb_mode, dataset, datum_path;
    bool        if_blocking, if_dualquant;
    double      mantissa, exponent;
    size_t      n_dim, d0, d1, d2, d3 = 1;

#if defined(_1D)
    cout << "\e[46mThis program is working for 1D datasets.\e[0m" << endl;
#elif defined(_2D)
    cout << "\e[46mThis program is working for 2D datasets.\e[0m" << endl;
#elif defined(_3D)
    cout << "\e[46mThis program is working for 3D datasets.\e[0m" << endl;
#endif

    if (argc != 8 && argc != 11) {
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
		if (argc == 8) {
            dataset = std::string(argv[6]);
            datum_path   = std::string(argv[7]);
        }
        else {
            n_dim        = atoi(&argv[6][1]);
            d0           = atoi(argv[7]);
            d1           = atoi(argv[8]);
            d2           = atoi(argv[9]);
            datum_path   = std::string(argv[10]);
        }

    }

    for_each(argv, argv + 8, [](auto i) { cout << i << " "; });
    cout << endl;
    auto eb_config  = new config_t(DICT_SIZE, mantissa, exponent);
    const int blk = BLK;
    auto dims_L16 = (argc == 8) ? InitializeDemoDims(dataset, DICT_SIZE, blk):InitializeDims(DICT_SIZE,n_dim,d0,d1,d2,d3,blk);
    printf("%-20s%s\n", "filename", datum_path.c_str());
    printf("%-20s%lu\n", "filesize", dims_L16[LEN] * sizeof(float));
    if (eb_mode == "r2r") {  // as of C++ 14, string is directly comparable?
        double value_range = GetDatumValueRange<float>(datum_path, dims_L16[LEN]);
        eb_config->ChangeToRelativeMode(value_range);
    }
    eb_config->debug();
    //    size_t c_byteSize;
    size_t num_outlier = 0;  // for calculating compression ratio

    cout << "block size:\t" << BLK << endl;
    auto ebs_L4 = InitializeErrorBoundFamily(eb_config);
    fm::cx_sim<float, int, BLK>(datum_path, dims_L16, ebs_L4, num_outlier, if_dualquant, if_blocking, true);
}
