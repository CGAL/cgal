// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//                 Matthias Baesken <baesken@informatik.uni-halle.de>

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#if defined (_MSC_VER)
#  include <sys/timeb.h>
#  include <sys/types.h>

#elif defined (__MINGW32__)
#  include <sys/timeb.h>
#  include <sys/types.h>

#else
// If none of the above PC compilers, use POSIX fct. gettimeofday()
#  include <sys/time.h>
#endif

namespace CGAL {

// Member functions for Real_timer
// =====================================

CGAL_INLINE_FUNCTION
double Real_timer::get_real_time() const {
    // Depends on the operating system.
    // Returns a (weakly ;-) monotone increasing time in seconds (with
    // possible wrap-around in case of overflow, see max()), or 0.0
    // if the system call for the time failed. If the system call
    // failed the static flag 'm_failed' is set and can be used
    // by the caller.
#if   defined(_MSC_VER)
    struct _timeb  t;
    _ftime_s(&t);
    return double(t.time) + double(t.millitm) / 1000.0;
#elif defined (__MINGW32__)
    struct timeb t;
    ftime(&t);
    return double(t.time) + double(t.millitm) / 1000.0;
#else // ! _MSC_VER && ! __MINGW32__//
    struct timeval t;
    int ret = gettimeofday( &t, nullptr);
    CGAL_warning_msg( ret == 0, "Call to gettimeofday() in class "
                      "CGAL::Real_timer failed - timings will be 0.");
    if ( ret == 0) {
        return double(t.tv_sec) + double(t.tv_usec) / 1000000;
    }
    get_static_realtimer_m_failed() = true;
    return 0.0;
#endif // ! _MSC_VER && ! __MINGW32__//
}

CGAL_INLINE_FUNCTION
double Real_timer::compute_precision() const {
    // Computes timer precision in seconds dynamically. Note that
    // the timer system call is probably non-trivial and will show
    // up in this time here (probably for one call). But that is just
    // fine that the call to the timer itself if reported as noise
    // in the precision.
    double min_res = DBL_MAX;
    for ( int i = 0; i < 5; ++i) {
        double current = get_real_time();
        if ( get_static_realtimer_m_failed() )
            return -1.0;
        double next    = get_real_time();
        while ( current >= next) { // wait until timer increases
            next = get_real_time();
            if ( get_static_realtimer_m_failed() )
                return -1.0;
        }
        // Get the minimum timing difference of all runs.
        if ( min_res > next - current)
            min_res = next - current;
    }
    return min_res;
}

CGAL_INLINE_FUNCTION
double Real_timer::precision() const {
    // computes precision upon first call
    // returns -1.0 if timer system call fails.
    static double prec = compute_precision();
    return prec;
}

} //namespace CGAL
