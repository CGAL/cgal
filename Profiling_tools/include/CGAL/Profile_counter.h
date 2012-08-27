// Copyright (c) 2005,2006,2008  INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_PROFILE_COUNTER_H
#define CGAL_PROFILE_COUNTER_H

// This file contains several classes to help in profiling, together with macros
// triggered by CGAL_PROFILE to enable them:
//
// - Profile_counter which is able to keep track of a number, and prints a
//   message in the destructor.  Typically, it can be used as a profile counter
//   in a static variable.
//
// - Profile_histogram_counter which is similar, but the counter is indexed by
//   a value (unsigned int), and the final dump is the histogram of the non-zero
//   counters.
//
// - Profile_branch_counter which keeps track of 2 counters, aiming at measuring
//   the ratio corresponding to the number of times a branch is taken.
//
// - Profile_branch_counter_3 which keeps track of 3 counters, aiming at measuring
//   the ratios corresponding to the number of times 2 branches are taken.
//
// See also CGAL/Profile_timer.h

// TODO :
// - Really complete the documentation!
// - Probably at some point we will need ways to selectively enable/disable profilers?
//   (kind-wise and/or place-wise)
// - Ideas for new kinds of profilers:
//   - lock counters in parallel mode
//     (e.g. time spent spinning, and/or number of locks taken or forbidden...)

#include <CGAL/config.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>

namespace CGAL {

struct Profile_counter
{
    Profile_counter(const std::string & ss)
      : i(0), s(ss) {}

    void operator++() { ++i; }

    ~Profile_counter()
    {
        std::cerr << "[CGAL::Profile_counter] "
                  << std::setw(10) << i << " " << s << std::endl;
    }

private:
    unsigned int i;
    const std::string s;
};



struct Profile_histogram_counter
{
    Profile_histogram_counter(const std::string & ss)
      : s(ss) {}

    void operator()(unsigned i) { ++counters[i]; }

    ~Profile_histogram_counter()
    {
        unsigned total=0;
        for (Counters::const_iterator it=counters.begin(), end=counters.end();
             it != end; ++it) {
            std::cerr << "[CGAL::Profile_histogram_counter] " << s;
            std::cerr << " [ " << std::setw(10) << it->first << " : "
                               << std::setw(10) << it->second << " ]"
                               << std::endl;
            total += it->second;
        }
        std::cerr << "[CGAL::Profile_histogram_counter] " << s;
        std::cerr << " [ " << std::setw(10) << "Total" << " : "
                           << std::setw(10) << total << " ]" << std::endl;
    }

private:
    typedef std::map<unsigned, unsigned>  Counters;
    Counters  counters;
    const std::string s;
};


struct Profile_branch_counter
{
    Profile_branch_counter(const std::string & ss)
      : i(0), j(0), s(ss) {}

    void operator++() { ++i; }

    void increment_branch() { ++j; }

    ~Profile_branch_counter()
    {
        std::cerr << "[CGAL::Profile_branch_counter] "
                  << std::setw(10) << j << " / "
                  << std::setw(10) << i << " " << s << std::endl;
    }

private:
    unsigned int i, j;
    const std::string s;
};


struct Profile_branch_counter_3
{
    Profile_branch_counter_3(const std::string & ss)
      : i(0), j(0), k(0), s(ss) {}

    void operator++() { ++i; }

    void increment_branch_1() { ++j; }
    void increment_branch_2() { ++k; }

    ~Profile_branch_counter_3()
    {
        std::cerr << "[CGAL::Profile_branch_counter_3] "
                  << std::setw(10) << k << " / "
                  << std::setw(10) << j << " / "
                  << std::setw(10) << i << " " << s << std::endl;
    }

private:
    unsigned int i, j, k;
    const std::string s;
};


#ifdef CGAL_PROFILE
#  define CGAL_PROFILER(Y) \
          { static CGAL::Profile_counter tmp(Y); ++tmp; }
#  define CGAL_HISTOGRAM_PROFILER(Y, Z) \
          { static CGAL::Profile_histogram_counter tmp(Y); tmp(Z); }
#  define CGAL_BRANCH_PROFILER(Y, NAME) \
          static CGAL::Profile_branch_counter NAME(Y); ++NAME;
#  define CGAL_BRANCH_PROFILER_BRANCH(NAME) \
          NAME.increment_branch();
#  define CGAL_BRANCH_PROFILER_3(Y, NAME) \
          static CGAL::Profile_branch_counter_3 NAME(Y); ++NAME;
#  define CGAL_BRANCH_PROFILER_BRANCH_1(NAME) \
          NAME.increment_branch_1();
#  define CGAL_BRANCH_PROFILER_BRANCH_2(NAME) \
          NAME.increment_branch_2();
#else
#  define CGAL_PROFILER(Y)
#  define CGAL_HISTOGRAM_PROFILER(Y, Z)
#  define CGAL_BRANCH_PROFILER(Y, NAME)
#  define CGAL_BRANCH_PROFILER_BRANCH(NAME)
#  define CGAL_BRANCH_PROFILER_3(Y, NAME) 
#  define CGAL_BRANCH_PROFILER_BRANCH_1(NAME)
#  define CGAL_BRANCH_PROFILER_BRANCH_2(NAME) 
#endif

} //namespace CGAL

#endif // CGAL_PROFILE_COUNTER_H
