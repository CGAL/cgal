// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>  
//                 Matthias Baesken <baesken@informatik.uni-halle.de>

#ifndef CGAL_TIMER_H
#define CGAL_TIMER_H 1

#include <CGAL/config.h>
#include <CGAL/assertions.h>

namespace CGAL {

// SECTION: A Timer for User-Process Time
// ========================================================================
// 
// DEFINITION
// 
// A timer `t' of type Timer is an object with a state. It is either
// running or it is stopped. The state is controlled with `t.start()'
// and `t.stop()'. The timer counts the time elapsed since its creation
// or last reset. It counts only the time where it is in the running
// state. The time information is given in seconds.

class CGAL_EXPORT Timer {
private:
    double      elapsed;
    double      started;
    int         interv;
    bool        running;

#ifdef CGAL_HEADER_ONLY
    static bool& get_static_timer_m_failed()
    {
      static bool m_failed = false;
      return m_failed;
    }
#else // CGAL_HEADER_ONLY
    static bool m_failed;
    static bool& get_static_timer_m_failed()
    { return CGAL::Timer::m_failed; }
#endif // CGAL_HEADER_ONLY

    double   user_process_time() const; // in seconds
    double   compute_precision() const; // in seconds
public:
    Timer() : elapsed(0.0), started(0.0), interv(0), running(false) {}

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
    double   max BOOST_PREVENT_MACRO_SUBSTITUTION ()        const;
};


// -----------------------------------------------------------------------

// Member functions for Timer
// ===========================

inline void Timer::start() {
    CGAL_precondition( ! running);
    started = user_process_time();
    running = true;
    ++ interv;
}

inline void Timer::stop() {
    CGAL_precondition( running);
    double t = user_process_time();
    elapsed += (t - started);
    started = 0.0;
    running = false;
}

inline void Timer::reset() {
    interv  = 0;
    elapsed = 0.0;
    if (running) {
	started = user_process_time();
	++ interv;
    } else {
        started = 0.0;
    }
}

inline double Timer::time() const {
    if (running) {
        double t = user_process_time();
	return elapsed + (t - started);
    }
    return elapsed;
}

} //namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Timer_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_TIMER_H //
// EOF //
