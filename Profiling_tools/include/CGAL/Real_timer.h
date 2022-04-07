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

#ifndef CGAL_REAL_TIMER_H
#define CGAL_REAL_TIMER_H 1

#include <CGAL/config.h>
#include <CGAL/assertions.h>
// For the numerical limits
#include <cfloat>

namespace CGAL {

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

class CGAL_EXPORT Real_timer {
private:
    double      elapsed;
    double      started;
    int         interv;
    bool        running;

#ifdef CGAL_HEADER_ONLY
    static bool& get_static_realtimer_m_failed()
    {
      static bool m_failed = false;
      return m_failed;
    }
#else // CGAL_HEADER_ONLY
    static bool m_failed;
    static bool& get_static_realtimer_m_failed()
    { return Real_timer::m_failed; }
#endif // CGAL_HEADER_ONLY

    double   get_real_time()     const; // in seconds
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
    double   max BOOST_PREVENT_MACRO_SUBSTITUTION ()        const { return DBL_MAX; }
};


// -----------------------------------------------------------------------

// Member functions for Real_timer
// ===========================

inline void Real_timer::start() {
    CGAL_precondition( ! running);
    started = get_real_time();
    running = true;
    ++ interv;
}

inline void Real_timer::stop() {
    CGAL_precondition( running);
    double t = get_real_time();
    elapsed += (t - started);
    started = 0.0;
    running = false;
}

inline void Real_timer::reset() {
    interv  = 0;
    elapsed = 0.0;
    if (running) {
        started = get_real_time();
        ++ interv;
    } else {
        started = 0.0;
    }
}

inline double Real_timer::time() const {
    if (running) {
        double t = get_real_time();
        return elapsed + (t - started);
    }
    return elapsed;
}

} //namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Real_timer_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_REAL_TIMER_H //
// EOF //
