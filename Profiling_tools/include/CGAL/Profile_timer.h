// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_PROFILE_TIMER_H
#define CGAL_PROFILE_TIMER_H

// This file contains a class to help in profiling, together with macros
// triggered by CGAL_PROFILE to enable it:
//
// - Profile_timer which keeps track of the time spent in its scope-block.
//   (be careful at recursive functions :-).
//
// See also CGAL/Profile_counter.h

#include <CGAL/config.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <CGAL/Real_timer.h>

CGAL_BEGIN_NAMESPACE

struct Profile_timer
{
    class Local {
        Profile_timer &p;
    public:
        Local(Profile_timer& p_) : p(p_) { p.start(); }
        ~Local() { p.stop(); }
    };

    Profile_timer(const std::string & ss)
      : s(ss) { t.reset(); }

    void start() { t.start(); }
    void stop() { t.stop(); }

    ~Profile_timer()
    {
        std::cerr << "[CGAL::Profile_timer] "
                  << std::setw(10) << t.time() << " " << s << std::endl;
    }

private:
    Real_timer t;
    const std::string s;
};


#ifdef CGAL_PROFILE
#  define CGAL_TIME_PROFILER(NAME) \
          static CGAL::Profile_timer CGAL_profile_timer_tmp(NAME); \
          CGAL::Profile_timer::Local CGAL_local_profile_timer_tmp(CGAL_profile_timer_tmp);
#else
#  define CGAL_TIME_PROFILER(NAME)
#endif

CGAL_END_NAMESPACE

#endif // CGAL_PROFILE_TIMER_H
