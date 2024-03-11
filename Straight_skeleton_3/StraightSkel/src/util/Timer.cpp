/**
 * @file   util/Timer.cpp
 * @author Gernot Walzl
 * @date   2013-06-25
 */

#include "util/Timer.h"

namespace util {

Timer::~Timer() {
    // intentionally does nothing
}

double Timer::now() {
    double result = 0.0;
    boost::posix_time::ptime now = boost::posix_time::microsec_clock::universal_time();
    boost::posix_time::time_duration duration = now - boost::posix_time::from_time_t(0);
    result = ((double)duration.total_milliseconds() / 1000.0);
    return result;
}

}
