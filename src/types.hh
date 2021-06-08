//
//  types.hh
//  waveSZ
//
//  Created by JianNan Tian on 6/8/19.
//  Copyright © 2019 JianNan Tian. All rights reserved.
//

#ifndef TYPES_HH
#define TYPES_HH

#include <algorithm>
#include <cmath>    // for FP32 bit representation
#include <cstddef>  // size_t
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "utils/format.hh"
#include "utils/io.hh"
#include "utils/timer.hh"

using namespace std;

template <typename T>
double GetDatumValueRange(string fname, size_t l);

void SetDims(size_t* dims_L16, size_t new_dims[4]);

typedef struct ErrorBoundConfigurator {
    int         capacity, radius;
    double      base, exp_base2, exp_base10;
    double      eb_base2, eb_base10, eb_final;
    std::string mode;

    void ChangeToRelativeMode(double value_range);

    void ChangeToTightBase2();

    ErrorBoundConfigurator(int _capacity = 32768, double _precision = 1, double eb = 0.0001);

    void debug() const;

} config_t;

// typedef struct DimensionInfo          dim_t;
// typedef struct ErrorBoundConfigurator config_t;

double* InitializeErrorBoundFamily(struct ErrorBoundConfigurator* eb_config);

typedef struct Stat {
    double minimum{}, maximum{}, range{};
    double PSNR{}, MSE{}, NRMSE{};
    double coeff{};
    double user_set_eb{}, max_abserr_vs_range{}, max_pwr_rel_abserr{};

    size_t len{}, max_abserr_index{};
    double max_abserr{};

} stat_t;

typedef struct Integer1  { int _0; }              Integer1;
typedef struct Integer2  { int _0, _1; }          Integer2;
typedef struct Integer3  { int _0, _1, _2; }      Integer3;
typedef struct Integer4  { int _0, _1, _2, _3; }  Integer4;
typedef struct UInteger1 { int _0; }             UInteger1;
typedef struct UInteger2 { int _0, _1; }         UInteger2;
typedef struct UInteger3 { int _0, _1, _2; }     UInteger3;
typedef struct UInteger4 { int _0, _1, _2, _3; } UInteger4;

#endif /* TYPES_HH */
