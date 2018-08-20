// Copyright (c) 2018 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_THREAD_H
#define CGAL_THREAD_H

#include <CGAL/config.h>

/*
  This file defines the following:
  - CGAL::cpp11::thread
  - CGAL::cpp11::atomic
  - CGAL::cpp11::sleep_for

  It uses either TBB or STD depending on what's available: as TBB can
  quite often override `std::thread`, it is possible that TBB will be
  used instead of STD even if the real CXX11 `std::thread` is
  available.

  As the conflicting API between TBB and STD can be quite complicated,
  we offer a more generic `sleep_for()` function that takes
  double-typed seconds as argument and deals with it.
*/

#if defined(CGAL_LINKED_WITH_TBB)
#  include <tbb/tbb_config.h>
#  if TBB_IMPLEMENT_CPP0X
#    include <tbb/compat/thread>
#    include <tbb/atomic.h>
#    include <tbb/tick_count.h>
#    define CGAL_USE_TBB_THREADS 1
#  else
#    define CGAL_USE_TBB_THREADS 0
#  endif
#else
#  define CGAL_USE_TBB_THREADS 0
#endif

#if !CGAL_USE_TBB_THREADS
#  include <thread>
#  include <chrono>
#endif

#include <CGAL/atomic.h> // for CGAL::cpp11::atomic

namespace CGAL {

namespace cpp11 {

#if CGAL_USE_TBB_THREADS
  
  using std::thread; // std::thread is declared by TBB if TBB_IMPLEMENT_CPP0X == 1

  inline void sleep_for (double seconds)
  {
    // std::this_thread::sleep_for is declared by TBB if TBB_IMPLEMENT_CPP0X == 1
    // It takes interval_t types as argument (!= from the std norm)
    std::this_thread::sleep_for(tbb::tick_count::interval_t(seconds));
  }
  
#else // C++11 implementation

  using std::thread;

  inline void sleep_for (double seconds)
  {
    // MSVC2013 cannot call `sleep_for()` with other types than
    // std::chrono::nanoseconds (bug in the implementation of the norm).
    typedef std::chrono::nanoseconds nanoseconds;
    nanoseconds ns (nanoseconds::rep (1000000000.0 * seconds));
    std::this_thread::sleep_for(ns);
  }

#endif

#if defined(CGAL_NO_ATOMIC) && defined(CGAL_LINKED_WITH_TBB)
  // If <CGAL/atomic.h> did not defined CGAL::cpp11::atomic, then use
  // tbb::atomic as a fallback.
  using tbb::atomic;
#endif

} // cpp11


} //namespace CGAL

#undef CGAL_USE_TBB_THREADS

#endif // CGAL_THREAD_H
