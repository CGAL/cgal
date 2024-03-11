/**
 * @file   util/Timer.h
 * @author Gernot Walzl
 * @date   2013-06-25
 */

#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

#include <boost/date_time/posix_time/posix_time.hpp>

namespace util {

class Timer {
public:
    virtual ~Timer();
    static double now();
protected:
    Timer();
};

}

#endif /* UTIL_TIMER_H */

