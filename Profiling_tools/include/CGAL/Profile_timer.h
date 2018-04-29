// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_PROFILE_TIMER_H
#define CGAL_PROFILE_TIMER_H

// This file contains a class to help in profiling, together with macros
// triggered by CGAL_PROFILE to enable it:
//
// - Profile_timer which keeps track of the time spent in its scope-block,
//   it accumulates the timings if it goes there several times.
//   Note: be careful at recursive functions...
//
// See also CGAL/Profile_counter.h
//
// TODO :
// - detect and handle recursion nicely (currently it accumulates recursive times).
// - parallel versions : per-thread and summing on all threads.
// - Timer vs Real_timer ?

#include <CGAL/config.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <CGAL/Timer.h>

#include <CGAL/disable_warnings.h>

namespace CGAL {

struct Profile_timer
{
    class Local {
        Timer rt;
        Profile_timer *p;
    public:
        Local(Profile_timer* p_) : p(p_) { rt.start(); }
        ~Local() { rt.stop(); p->t += rt.time(); }
    };

    Profile_timer(const std::string & ss)
      : t(0), s(ss) {}

    ~Profile_timer()
    {
        std::cerr << "[CGAL::Profile_timer] "
                  << std::setw(10) << t << " " << s << std::endl;
    }

private:
    double t;
    const std::string s;
};


#ifdef CGAL_PROFILE
#  define CGAL_TIME_PROFILER(NAME) \
          static CGAL::Profile_timer CGAL_profile_timer_tmp(NAME); \
          CGAL::Profile_timer::Local CGAL_local_profile_timer_tmp(&CGAL_profile_timer_tmp);
#else
#  define CGAL_TIME_PROFILER(NAME)
#endif

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_PROFILE_TIMER_H
