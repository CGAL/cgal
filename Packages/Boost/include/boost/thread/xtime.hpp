// Copyright (C) 2001-2003
// William E. Kempf
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  William E. Kempf makes no representations
// about the suitability of this software for any purpose.
// It is provided "as is" without express or implied warranty.

#ifndef BOOST_XTIME_WEK070601_HPP
#define BOOST_XTIME_WEK070601_HPP

#include <boost/thread/detail/config.hpp>

#include <boost/cstdint.hpp>

namespace boost {

enum xtime_clock_types
{
    TIME_UTC=1
//    TIME_TAI,
//    TIME_MONOTONIC,
//    TIME_PROCESS,
//    TIME_THREAD,
//    TIME_LOCAL,
//    TIME_SYNC,
//    TIME_RESOLUTION
};

struct xtime
{
#if defined(BOOST_NO_INT64_T)
    typedef int_fast32_t xtime_sec_t; //INT_FAST32_MIN <= sec <= INT_FAST32_MAX
#else
    typedef int_fast64_t xtime_sec_t; //INT_FAST64_MIN <= sec <= INT_FAST64_MAX
#endif

    typedef int_fast32_t xtime_nsec_t; //0 <= xtime.nsec < NANOSECONDS_PER_SECOND

    xtime_sec_t sec;
    xtime_nsec_t nsec;
};

int BOOST_THREAD_DECL xtime_get(struct xtime* xtp, int clock_type);

inline int xtime_cmp(const xtime& xt1, const xtime& xt2)
{
    if (xt1.sec == xt2.sec)
        return (int)(xt1.nsec - xt2.nsec);
    else 
        return (xt1.sec > xt2.sec) ? 1 : -1;
}

} // namespace boost

#endif //BOOST_XTIME_WEK070601_HPP
