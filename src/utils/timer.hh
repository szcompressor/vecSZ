//
// Created by JianNan Tian on 2019-08-26.
//

#ifndef TIMER_HH
#define TIMER_HH

#include <chrono>

using std::cerr;
using std::cout;
using std::endl;

using hires = std::chrono::high_resolution_clock;
typedef std::chrono::duration<double>                               duration_t;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> hires_clock_t;

/**
 * @brief A timer wrapper for arbitrary function (no handling return for now);
 * Adapted from https://stackoverflow.com/a/33900479/8740097 (CC BY-SA 3.0)
 *
 * @tparam F auto function type
 * @tparam Args variadic function argument type
 * @param func non-return function to be timed
 * @param args variadic function arguments
 * @return double time in seconds
 */
template <typename F, typename... Args>
double TimeThisFunction(F func, Args&&... args)
{
    auto t0 = hires::now();
    func(std::forward<Args>(args)...);
    return static_cast<duration_t>(hires::now() - t0).count();
}

#endif  // TIMER_HH
