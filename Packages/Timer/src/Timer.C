// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : src/Timer.C
// package       : Timer (1.9)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.9
// revision_date : 4 April 2001 
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>  
//                 Matthias Baesken <baesken@informatik.uni-halle.de>
// coordinator   : INRIA, Sophia Antipolis
//
// A timer class to measure cpu time of user-processes.
// ======================================================================

#include <CGAL/Timer.h>

// Determine if the POSIX function getrusage is available, otherwise
// use the previous solution based on std::clock().
// First, detect POSIX. We cannot reliably use "unistd.h", 
// but limits.h is part of the C standard.
#include <climits>
#ifdef _POSIX_ARG_MAX // now that should be POSIX
#include <unistd.h>
#ifdef _POSIX_VERSION
#ifdef _XOPEN_UNIX // XSI: X/Open System Interfaces Extension
#define CGAL__GETRUSAGE 1
#endif
#endif
#endif


#ifdef CGAL__GETRUSAGE
// types, function prototype and constants for the POSIX function
// int getrusage (int who, struct rusage *usage);
#include <sys/resource.h>
// For the numerical limits
#else //  CGAL__GETRUSAGE //
// used for clock()
#include <ctime>
#endif //  CGAL__GETRUSAGE //

// For the numerical limits
#include <float.h>


CGAL_BEGIN_NAMESPACE


// Static member variable for Timer
// =====================================

bool Timer::m_failed = false;

// Member functions for Timer
// =====================================

double Timer::user_process_time() const {
    // Depends on the operating system.
    // Returns a (weakly ;-) monotone increasing time in seconds (with
    // possible wrap-around in case of overflow, see max()), or 0.0
    // if the system call for the time failed. If the system call
    // failed the static flag 'm_failed' is set and can be used
    // by the caller.
#ifdef CGAL__GETRUSAGE
    struct rusage usage;
    int ret = getrusage( RUSAGE_SELF, &usage);
    CGAL_warning_msg( ret == 0, "Call to getrusage() in class CGAL::Timer "
                      "failed - timings will be 0.");
    if ( ret == 0) {
        return double( usage.ru_utime.tv_sec)               // seconds
             + double( usage.ru_utime.tv_usec) / 1000000.0; // microseconds
    }
#else // CGAL__GETRUSAGE //
    CGAL_CLIB_STD::clock_t clk = CGAL_CLIB_STD::clock();
    CGAL_warning_msg( clk != (CGAL_CLIB_STD::clock_t)-1,
        "Call to clock() in class CGAL::Timer failed - timings will be 0.");
    if ( clk != (CGAL_CLIB_STD::clock_t)-1) {
        return double(clk) / CLOCKS_PER_SEC;
    }        
#endif // CGAL__GETRUSAGE //
    m_failed = true;
    return 0.0;
}

double Timer::compute_precision() const {
    // Computes timer precision in seconds dynamically. Note that
    // the timer system call is probably non-trivial and will show
    // up in this time here (probably for one call). But that is just
    // fine that the call to the timer itself if reported as noise 
    // in the precision.
    double min_res = DBL_MAX;
    for ( int i = 0; i < 5; ++i) {
        double current = user_process_time();
        if ( m_failed)
            return -1.0;
        double next    = user_process_time();
        while ( current >= next) { // wait until timer increases
            next = user_process_time();
            if ( m_failed)
                return -1.0;
        }
        // Get the minimum timing difference of all runs.
        if ( min_res > next - current)
            min_res = next - current;
    }
    return min_res;
}

double Timer::precision() const {
    // computes precision upon first call
    // returns -1.0 if timer system call fails.
    static double prec = compute_precision();
    return prec;
}

double Timer::max() const { 
    // Depends on the operating system.
#ifdef CGAL__GETRUSAGE
    return DBL_MAX;
#else // CGAL__GETRUSAGE //
    return 2146.0;
#endif // CGAL__GETRUSAGE //
}

CGAL_END_NAMESPACE

// EOF //
