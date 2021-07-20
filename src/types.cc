//
//  types.hh
//  waveSZ
//
//  Created by JianNan Tian on 6/8/19.
//  Copyright © 2019 JianNan Tian. All rights reserved.
//

#include <algorithm>
#include <cmath>    // for FP32 bit representation
#include <cstddef>  // size_t
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "constants.hh"
#include "utils/format.hh"
#include "utils/io.hh"
#include "utils/timer.hh"
#include "types.hh"

using namespace std;

template <typename T>
double GetDatumValueRange(string fname, size_t l)
{
    auto d    = io::ReadBinaryToNewArray<T>(fname, l);
    T    max_ = *std::max_element(d, d + l);
    T    min_ = *std::min_element(d, d + l);
    delete[] d;
    return max_ - min_;
}

template double GetDatumValueRange<float>(string fname, size_t l);
template double GetDatumValueRange<double>(string fname, size_t l);

ErrorBoundConfigurator::ErrorBoundConfigurator(int _capacity, double eb)
{
    int _base = 10;
    capacity = _capacity;
    radius   = capacity / 2;
    mode     = std::string("ABS");

    eb_final   = eb;
}

void ErrorBoundConfigurator::ChangeToRelativeMode(double value_range)
{
    if (value_range == 0) {
        cerr << log_err << "INVALID VALUE RANGE!" << endl;
        exit(1);
    }
    cout << log_info << "change to r2r mode \e[2m(relative-to-value-range)\e[0m" << endl;
    cout << log_null << "eb --> " << eb_final << " x " << value_range << " = ";
    this->eb_final *= value_range;
    cout << eb_final << endl;
    mode = std::string("VRREL");
}

void ErrorBoundConfigurator::ChangeToTightBase2()
{
    base = 2;
    cout << log_info << "switch.to.tight.base2.mode, eb changed from " << eb_final << " = 2^(" << exp_base2 << ") to ";
    cout << "the exp base2 before changing:\t" << exp_base2 << endl;
    exp_base2 = floor(exp_base2);
    cout << "the exp base2 after changing:\t" << exp_base2 << endl;
    eb_final = pow(2, exp_base2);
    cout << eb_final << " = 2^(" << exp_base2 << ")" << endl;
}

void ErrorBoundConfigurator::debug() const
{
    cout << log_dbg << "exp.base10:\t" << exp_base10 << endl;
    cout << log_dbg << "exp.base2:\t" << exp_base2 << endl;
    cout << log_dbg << "final.eb:\t" << eb_final << endl;
}

//} config_t;

typedef struct ErrorBoundConfigurator config_t;

double* InitializeErrorBoundFamily(config_t* eb_config)
{
    auto ebs_L4 = new double[4]();
    ebs_L4[0]   = eb_config->eb_final;            // eb
    ebs_L4[1]   = 1 / eb_config->eb_final;        // 1/eb
    ebs_L4[2]   = 2 * eb_config->eb_final;        // 2* eb
    ebs_L4[3]   = 1 / (2 * eb_config->eb_final);  // 1/(2*eb)
    cout << log_dbg << "final.eb:\t" << eb_config->eb_final << endl;
    return ebs_L4;
}
