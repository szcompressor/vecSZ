// 200209

#ifndef ANALYSIS_HH
#define ANALYSIS_HH

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "../constants.hh"
#include "format.hh"

#define TOTAL_BORDER 0
#define TOP_BORDER 1
#define LEFT_BORDER 2
#define FRONT_BORDER 3

using std::cerr;
using std::cout;
using std::endl;
using std::tuple;

namespace analysis {
template <typename T>
tuple<double, double, double> getStat(T* __d, size_t l, bool print = false)
{
    double Min = *std::min_element(__d, __d + l);
    double Max = *std::max_element(__d, __d + l);
    double sum = std::accumulate(__d, __d + l, 0);
    double rng = Max - Min;
    double avg = sum / l;
    if (print) {
        cout << "rng: " << rng << endl;
        cout << "min: " << Min << endl;
        cout << "max: " << Max << endl;
        cout << "avg: " << avg << endl;
    }
    return std::make_tuple(Max, Min, rng);
}

template <typename T>
tuple<double, double, double> getStat(std::vector<T> __d, bool print = false)
{
    double Min = *std::min_element(__d.begin(), __d.end());
    double Max = *std::max_element(__d.begin(), __d.end());
    double sum = std::accumulate(__d.begin(), __d.end(), 0);
    double rng = Max - Min;
    double avg = sum / __d.size();
    if (print) {
        cout << "rng: " << rng << endl;
        cout << "min: " << Min << endl;
        cout << "max: " << Max << endl;
        cout << "avg: " << avg << endl;
    }
    return std::make_tuple(Max, Min, rng);
}

template <typename T>
void getEntropy(T* code, size_t l, size_t cap = 1024)
{
    if (cap == 0) {
        cerr << "wrong cap" << endl;
        exit(1);
    }
    auto arr = new size_t[cap]();
    for (size_t i = 0; i < l; i++) arr[code[i]]++;
    std::vector<double> raw(arr, arr + cap);
    std::vector<double> frequencies;
    std::copy_if(raw.begin(), raw.end(), std::back_inserter(frequencies), [](double& e) { return e != 0; });
    double entropy = 0;
    for (auto freq : frequencies) {
        //        cout << -(freq / l) * log2(freq / l) << endl;
        entropy += -(freq / l) * log2(freq / l);
    }

    cout << "entropy:\t" << entropy << endl;
    delete[] arr;
}

// TODO automatically omit bins that are less than 1%
template <typename T>
void histogram(
    const std::string& tag,
    T*                 __d_POD,
    size_t             l,
    size_t             __bins                  = 16,
    bool               log_freq                = false,
    double             override_min            = 0,
    double             override_max            = 0,
    bool               eliminate_zeros         = false,
    bool               use_scientific_notation = true)
{
    std::vector<T> __d(__d_POD, __d_POD + l);
    std::vector<T> __d_nonzero;
    //    std::vector<size_t> arr;
    //    arr.reserve(__bins);
    //    for (size_t i = 0; i< __bins; i++) arr.push_back(0);
    auto arr = new size_t[__bins]();

    if (eliminate_zeros) {
        std::copy_if(__d.begin(), __d.end(), std::back_inserter(__d_nonzero), [](int i) { return i != 0; });
    }
    double Min = *std::min_element(__d.begin(), __d.end());
    double Max = *std::max_element(__d.begin(), __d.end());
    //    double sum = std::accumulate(__d.begin(), __d.end(), 0);
    double rng = Max - Min;
    //    double avg = sum / l;

    cout << "\e[7m[[" << tag << "]]\e[0m";
    if (override_max > override_min) {
        cout << "zoom into " << override_min << "--" << override_max << endl;
        std::tie(Max, Min, rng) = std::make_tuple(override_max, override_min, override_max - override_min);
    }
    double step = rng / __bins;
    for (size_t i = 0; i < l; i++) arr[static_cast<size_t>((__d[i] - Min) / step)]++;
    std::vector<size_t> __viz(arr, arr + __bins);
    //    std::vector<size_t> __viz(arr);

    // visualization
    printf("\tbins:\t%zu\tbin_width:\t%lf\n", __bins, step);
    //    printf("count:\t%zu\tmin:\t%lf\tmax:\t%lf\trng:\t%lf\n", l, Min, Max, rng);
    cout << "count:\t" << l << "\t";
    cout << "min:\t" << Min << "\t";
    cout << "max:\t" << Max << "\t";
    cout << "rng:\t" << rng << endl;

    if (log_freq) {
        cout << "using log_freq" << endl;
        std::for_each(__viz.begin(), __viz.end(), [](size_t& n) { n = log2(n); });
    }

    size_t longest     = *std::max_element(__viz.begin(), __viz.end());
    size_t bar_str_len = 64;  // scale according to the longest
    std::for_each(__viz.begin(), __viz.end(), [&](size_t& n) { n = static_cast<size_t>(n / static_cast<double>(longest) * bar_str_len); });

    for (size_t i = 0; i < __bins; i++) {
        // normalize to width
        cout << "|"
             << "\33[43m";

        for (size_t j = 0; j < bar_str_len + 1; j++) {
            if (j < __viz[i]) cout << "-";
            else if (j == __viz[i])
                cout << "\33[0m"
                     << "+";
            else
                cout << " ";
        }
        cout.precision(2);
        cout << "    ";
        if (use_scientific_notation) cout << std::scientific;
        cout << Min + i * step << " -- " << Min + (i + 1) * step;
        cout << "  ";
        cout << std::setw((int)log10(l) + 2);
        cout << arr[i];
        cout << "   ";
        cout << std::defaultfloat << std::setw(5) << arr[i] / static_cast<double>(l) * 100 << "%" << endl;
    }
    cout << endl;
    //    delete[] arr;
}

template<typename C>
size_t* get_border_outliers_1d(size_t* outlier, C* code, size_t* dims, int blksz, int b0)
{
	size_t idx0 = b0 * blksz;
	for (size_t i0 = 0; i0 < blksz; i0++)
	{
		size_t id = idx0 + i0;
		if (id >= dims[DIM0]) continue;
		if (code[id] == 0)
		{
			outlier[LEFT_BORDER] += (i0 == 0) ? 1 : 0;
			outlier[TOTAL_BORDER] += (i0 == 0) ? 1 : 0;
		}

	}

	return outlier;
}

template<typename C>
size_t* get_border_outliers_2d(size_t* outlier, C* code, size_t* dims, int blksz, int b0, int b1)
{
	size_t idx0 = b0 * blksz;
	size_t idx1 = b1 * blksz;
	for (size_t i1 = 0; i1 < blksz; i1++)
	{
		for (size_t i0 = 0; i0 < blksz; i0++)
		{
			size_t gi1 = idx1 + i1;
			size_t gi0 = idx0 + i0;
			if (gi1 >= dims[DIM1] or gi0 >= dims[DIM0]) continue;
			size_t id  = gi0 + gi1 * dims[DIM0];
			if (code[id] == 0)
			{
				outlier[TOP_BORDER]  += (i0 == 0) ? 1 : 0;
				outlier[LEFT_BORDER] += (i1 == 0) ? 1 : 0;
				outlier[TOTAL_BORDER] += ((i0 == 0) || (i1 == 0)) ? 1 : 0;
			}
		}
	}

	return outlier;
}

template<typename C>
size_t* get_border_outliers_3d(size_t* outlier, C* code, size_t* dims, int blksz, int b0, int b1, int b2)
{
	size_t idx0 = b0 * blksz;
	size_t idx1 = b1 * blksz;
	size_t idx2 = b2 * blksz;
	for (size_t i2 = 0; i2 < blksz; i2++)
	{
		for (size_t i1 = 0; i1 < blksz; i1++)
		{
			for (size_t i0 = 0; i0 < blksz; i0++)
			{
				size_t gi2 = idx2 + i2;
				size_t gi1 = idx1 + i1;
				size_t gi0 = idx0 + i0;
				if (gi2 >= dims[DIM2] or gi1 >= dims[DIM1] or gi0 >= dims[DIM0]) continue;
				size_t id  = gi0 + gi1 * dims[DIM0] + gi2 * dims[DIM1] * dims[DIM0];
				if (code[id] == 0)
				{
					outlier[TOP_BORDER]  += (i0 == 0) ? 1 : 0;
					outlier[LEFT_BORDER] += (i1 == 0) ? 1 : 0;
					outlier[FRONT_BORDER] += (i2 == 0) ? 1 : 0;
					outlier[TOTAL_BORDER] += ((i0 == 0) || (i1 == 0) || (i2 == 0)) ? 1 : 0;
				}
			}
		}
	}

	return outlier;
}

template<typename C>
void get_outliers(size_t* dims, C* code, int blksz, size_t numOutlier)
{
    size_t* outlier = new size_t[4];
    outlier[TOTAL_BORDER] = 0;
    outlier[TOP_BORDER] = 0;
    outlier[LEFT_BORDER] = 0;
    outlier[FRONT_BORDER] = 0;

    if (dims[nDIM] == 1)
    {
        #pragma omp parallel for
        for (size_t b0 = 0; b0 < dims[nBLK0]; b0++)
        {
	    outlier = get_border_outliers_1d(outlier, code, dims, blksz, b0);
        }
    }
    else if (dims[nDIM] == 2)
    {
        #pragma omp parallel for
        for (size_t b1 = 0; b1 < dims[nBLK1]; b1++)
        {
            for (size_t b0 = 0; b0 < dims[nBLK0]; b0++)
            {
	    	outlier = get_border_outliers_2d(outlier, code, dims, blksz, b0, b1);
            }
        }
    }
    else if (dims[nDIM] == 3)
    {
        #pragma omp parallel for
        for (size_t b2 = 0; b2 < dims[nBLK2]; b2++)
        {
            for (size_t b1 = 0; b1 < dims[nBLK1]; b1++)
            {
                for (size_t b0 = 0; b0 < dims[nBLK0]; b0++)
                {
	    	    outlier = get_border_outliers_3d(outlier, code, dims, blksz, b0, b1, b2);
                }
            }
        }
    }

    float totBorder  = (outlier[TOTAL_BORDER] == 0) ? 0 : ( (float) outlier[TOTAL_BORDER] / numOutlier ) * 100;
    float leftBorder = (outlier[LEFT_BORDER]  == 0) ? 0 : ( (float) outlier[LEFT_BORDER]  / numOutlier ) * 100;
    float topBorder  = (outlier[TOP_BORDER]   == 0) ? 0 : ( (float) outlier[TOP_BORDER]   / numOutlier ) * 100;
    float frntBorder = (outlier[FRONT_BORDER] == 0) ? 0 : ( (float) outlier[FRONT_BORDER] / numOutlier ) * 100;

    // report results
    cout << "--------------------- Outliers ---------------------" << endl;
    cout << log_info << "Total Outliers:      " << numOutlier << "  " << endl;
    cout << log_info << "\% Outliers:         " << (float) numOutlier / (dims[LEN]) << " \%" << endl;
    cout << log_info << "Border Outliers:     " << outlier[TOTAL_BORDER] << "  " << endl;
    cout << log_info << "\% Border Outliers:  " << totBorder << " \%"  << endl;
    cout << log_info << "  Left Outliers:     " << outlier[LEFT_BORDER] << "  " << endl;
    cout << log_info << "  \% Left Outliers:  " << leftBorder << " \%"  << endl;
    cout << log_info << "  Top Outliers:      " << outlier[TOP_BORDER] << "  " << endl;
    cout << log_info << "  \% Top Outliers:   " << topBorder << " \%"  << endl;
    cout << log_info << "  Front Outliers:    " << outlier[FRONT_BORDER] << "  " << endl;
    cout << log_info << "  \% Front Outliers: " << frntBorder << " \%"  << endl;
    cout << "----------------------------------------------------" << endl;
}

}  // namespace analysis

#endif
