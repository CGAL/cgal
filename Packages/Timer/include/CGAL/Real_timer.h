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
// file          : include/CGAL/Real_timer.h
// package       : Timer (1.9)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.9
// revision_date : 4 April 2001 
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>  
//                 Matthias Baesken <baesken@informatik.uni-halle.de>
// coordinator   : INRIA, Sophia Antipolis
//
// A timer class to measure real-time.
// ======================================================================

#ifndef CGAL_REAL_TIMER_H
#define CGAL_REAL_TIMER_H 1

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

// SECTION: A Timer Measuring Real-Time
// ========================================================================
// 
// DEFINITION
// 
// A timer `t' of type Real_timer is an object with a state. It is either
// running or it is stopped. The state is controlled with `t.start()'
// and `t.stop()'. The timer counts the time elapsed since its creation
// or last reset. It counts only the time where it is in the running
// state. The time information is given in seconds.

class Real_timer {
private:
    double      elapsed;
    double      started;
    int         interv;
    bool        running;

    static bool m_failed;

    double   user_process_time() const; // in seconds
    double   compute_precision() const; // in seconds
public:
    Real_timer() : elapsed(0.0), started(0.0), interv(0), running(false) {}

    void     start();
    void     stop ();
    void     reset();
    bool     is_running() const { return running; }

    double   time()       const;
    int      intervals()  const { return interv; }
    double   precision()  const;
    // Returns timer precison. Computes it dynamically at first call.
    // Returns -1.0 if timer system call fails, which, for a proper coded
    // test towards precision leads to an immediate stop of an otherwise 
    // infinite loop (fixed tolerance * total time >= precision).
    double   max()        const { return DBL_MAX; }
};


// -----------------------------------------------------------------------

// Member functions for Real_timer
// ===========================

inline void Real_timer::start() {
    CGAL_precondition( ! running);
    started = user_process_time();
    running = true;
    ++ interv;
}

inline void Real_timer::stop() {
    CGAL_precondition( running);
    elapsed += user_process_time() - started;
    started = 0.0;
    running = false;
}

inline void Real_timer::reset() {
    interv  = 0;
    elapsed = 0.0;
    if (running) {
	started = user_process_time();
	++ interv;
    } else {
        started = 0.0;
    }
}

inline double Real_timer::time() const {
    if (running)
	return elapsed + (user_process_time() - started);
    return elapsed;
}

CGAL_END_NAMESPACE

#endif // CGAL_REAL_TIMER_H //
// EOF //
