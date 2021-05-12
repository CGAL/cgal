//r Copyright (c) 2005,2006,2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
//  If CGAL_CONCURRENT_PROFILE is defined, the counters can be concurrently updated
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
#include <sstream>
#include <iomanip>
#include <string>
#include <map>

#include <CGAL/disable_warnings.h>

// Automatically define CGAL_CONCURRENT_PROFILE if we're linked with TBB
#ifdef CGAL_LINKED_WITH_TBB
# ifndef CGAL_CONCURRENT_PROFILE
#   define CGAL_CONCURRENT_PROFILE
# endif
#else
// Automatically UNdefine CGAL_CONCURRENT_PROFILE if we're NOT linked with TBB
# ifdef CGAL_CONCURRENT_PROFILE
#   undef CGAL_CONCURRENT_PROFILE
# endif
#endif

#ifdef CGAL_CONCURRENT_PROFILE
# include "tbb/concurrent_hash_map.h"
#endif

namespace CGAL {

  namespace internal {
    struct dotted : std::numpunct<char> {
      char do_thousands_sep()   const { return '.'; }  // separate with dots
      std::string do_grouping() const { return "\3"; } // groups of 3 digits
      static void imbue(std::ostream &os) {
        os.imbue(std::locale(os.getloc(), new dotted));
      }
    };

    inline std::string dot_it(int i)
    {
      std::stringstream ss;
      dotted::imbue(ss);
      ss << i;
      return ss.str();
    }
  }

struct Profile_counter
{
    Profile_counter(const std::string & ss)
      : s(ss)
    {
      i = 0; // needed here because of tbb::atomic
    }

    void operator++() { ++i; }

    ~Profile_counter()
    {


      std::cerr << "[CGAL::Profile_counter] "
                << std::setw(10) << internal::dot_it(i) << " " << s << std::endl;
    }

private:
#ifdef CGAL_CONCURRENT_PROFILE
    tbb::atomic<unsigned int> i;
#else
    unsigned int i;
#endif
    const std::string s;
};



struct Profile_histogram_counter
{
private:
#ifdef CGAL_CONCURRENT_PROFILE
    typedef tbb::concurrent_hash_map<unsigned, unsigned>  Counters;
#else
    typedef std::map<unsigned, unsigned>  Counters;
#endif

public:
    Profile_histogram_counter(const std::string & ss)
      : s(ss) {}

    void operator()(unsigned i)
    {
#ifdef CGAL_CONCURRENT_PROFILE
      Counters::accessor a;
      counters.insert(a, i);
      ++a->second;
#else
      ++counters[i];
#endif
    }

    ~Profile_histogram_counter()
    {
        unsigned total=0;
        for (Counters::const_iterator it=counters.begin(), end=counters.end();
             it != end; ++it) {
            std::cerr << "[CGAL::Profile_histogram_counter] " << s;
            std::cerr << " [ " << std::setw(10) << internal::dot_it(it->first) << " : "
                      << std::setw(10) << internal::dot_it(it->second) << " ]"
                               << std::endl;
            total += it->second;
        }
        std::cerr << "[CGAL::Profile_histogram_counter] " << s;
        std::cerr << " [ " << std::setw(10) << "Total" << " : "
                           << std::setw(10) << total << " ]" << std::endl;
    }

private:
    Counters  counters;
    const std::string s;
};


struct Profile_branch_counter
{
    Profile_branch_counter(const std::string & ss)
      : s(ss)
    {
      i = j = 0; // needed here because of tbb::atomic
    }

    void operator++() { ++i; }

    void increment_branch() { ++j; }

    ~Profile_branch_counter()
    {
        std::cerr << "[CGAL::Profile_branch_counter] "
                  << std::setw(10) << internal::dot_it(j) << " / "
                  << std::setw(10) << internal::dot_it(i) << " " << s << std::endl;
    }

private:
#ifdef CGAL_CONCURRENT_PROFILE
    tbb::atomic<unsigned int> i, j;
#else
    unsigned int i, j;
#endif
    const std::string s;
};


struct Profile_branch_counter_3
{
    Profile_branch_counter_3(const std::string & ss)
      : s(ss)
    {
      i = j = k = 0; // needed here because of tbb::atomic
    }

    void operator++() { ++i; }

    void increment_branch_1() { ++j; }
    void increment_branch_2() { ++k; }

    ~Profile_branch_counter_3()
    {
        std::cerr << "[CGAL::Profile_branch_counter_3] "
                  << std::setw(10) << internal::dot_it(k) << " / "
                  << std::setw(10) << internal::dot_it(j) << " / "
                  << std::setw(10) << internal::dot_it(i) << " " << s << std::endl;
    }

private:
#ifdef CGAL_CONCURRENT_PROFILE
    tbb::atomic<unsigned int> i, j, k;
#else
    unsigned int i, j, k;
#endif
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

#include <CGAL/enable_warnings.h>

#endif // CGAL_PROFILE_COUNTER_H
