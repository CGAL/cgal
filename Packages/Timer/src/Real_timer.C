// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>  
//                 Matthias Baesken <baesken@informatik.uni-halle.de>

#include <CGAL/Real_timer.h>

#if defined (_MSC_VER)
#define CGAL_PROTECT_SYS_TIME_H
#include <sys/timeb.h>
#include <sys/types.h>
#endif

#if defined (__BORLANDC__)
#define CGAL_PROTECT_SYS_TIME_H
#include <sys/timeb>
#include <ctype> 
#endif

#if defined (__MWERKS__)
#define CGAL_PROTECT_SYS_TIME_H
#include <ctime>
#endif

#if defined (__MINGW32__)
#define CGAL_PROTECT_SYS_TIME_H
  #include <sys/timeb.h>
  #include <sys/types.h>
#endif 

// If none of the above PC compilers, use POSIX fct. gettimeofday()
#ifndef CGAL_PROTECT_SYS_TIME_H
#include <sys/time.h>
#define CGAL_PROTECT_SYS_TIME_H
#endif // CGAL_PROTECT_SYS_TIME_H


CGAL_BEGIN_NAMESPACE


// Static member variable for Real_timer
// =====================================

bool Real_timer::m_failed = false;

// Member functions for Real_timer
// =====================================

double Real_timer::get_real_time() const {
    // Depends on the operating system.
    // Returns a (weakly ;-) monotone increasing time in seconds (with
    // possible wrap-around in case of overflow, see max()), or 0.0
    // if the system call for the time failed. If the system call
    // failed the static flag 'm_failed' is set and can be used
    // by the caller.
#if   defined(_MSC_VER)
    struct _timeb  t;
    _ftime(&t);  
    return double(t.time) + double(t.millitm) / 1000.0;
#elif defined(__BORLANDC__)
    struct timeb t;
    ftime(&t);  
    return double(t.time) + double(t.millitm) / 1000.0;
#elif defined(__MWERKS__)
    static bool initialized = false;
    static CGAL_CLIB_STD::time_t t_begin;
    if ( ! initialized) {
        CGAL_CLIB_STD::time_t tt;
        CGAL_CLIB_STD::time(&tt);
        struct std::tm* lt = CGAL_CLIB_STD::localtime(&tt);
        /// we copy them because mktime can modify
        struct std::tm t = *lt;
        t_begin = CGAL_CLIB_STD::mktime(&t);
    }
    CGAL_CLIB_STD::time_t tt;
    CGAL_CLIB_STD::time(&tt);
    struct std::tm* lt = CGAL_CLIB_STD::localtime(&tt);
    /// we copy them because mktime can modify
    struct std::tm t = *lt;
    CGAL_CLIB_STD::time_t t_end = CGAL_CLIB_STD::mktime(&t);
    return CGAL_CLIB_STD::difftime(t_end, t_begin);
#elif defined (__MINGW32__)
    struct timeb t;
    ftime(&t);
    return double(t.time) + double(t.millitm) / 1000.0;
#else // ! _MSC_VER && ! __BORLANDC__ && ! __MWERKS__ && ! __MINGW32__//
    struct timeval t;
    int ret = gettimeofday( &t, NULL);
    CGAL_warning_msg( ret == 0, "Call to gettimeofday() in class "
                      "CGAL::Real_timer failed - timings will be 0.");
    if ( ret == 0) {
        return double(t.tv_sec) + double(t.tv_usec) / 1000000;
    }
    m_failed = true;
    return 0.0;
#endif // ! _MSC_VER && ! __BORLANDC__ && ! __MWERKS__  && ! __MINGW32__//
}

double Real_timer::compute_precision() const {
    // Computes timer precision in seconds dynamically. Note that
    // the timer system call is probably non-trivial and will show
    // up in this time here (probably for one call). But that is just
    // fine that the call to the timer itself if reported as noise 
    // in the precision.
    double min_res = DBL_MAX;
    for ( int i = 0; i < 5; ++i) {
        double current = get_real_time();
        if ( m_failed)
            return -1.0;
        double next    = get_real_time();
        while ( current >= next) { // wait until timer increases
            next = get_real_time();
            if ( m_failed)
                return -1.0;
        }
        // Get the minimum timing difference of all runs.
        if ( min_res > next - current)
            min_res = next - current;
    }
    return min_res;
}

double Real_timer::precision() const {
    // computes precision upon first call
    // returns -1.0 if timer system call fails.
    static double prec = compute_precision();
    return prec;
}

CGAL_END_NAMESPACE

// EOF //
