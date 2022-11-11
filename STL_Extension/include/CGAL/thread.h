// Copyright (c) 2018 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_THREAD_H
#define CGAL_THREAD_H

#define CGAL_DEPRECATED_HEADER "<CGAL/thread.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/config.h>"

#include <CGAL/config.h>

/*
  This file defines the following:
  - CGAL::cpp11::thread
  - CGAL::cpp11::atomic
  - CGAL::cpp11::sleep_for
*/

#include <thread>
#include <chrono>
#include <CGAL/atomic.h> // for CGAL::cpp11::atomic

namespace CGAL {

namespace cpp11 {

  using std::thread;

  inline void sleep_for (double seconds)
  {
    // MSVC2013 cannot call `sleep_for()` with other types than
    // std::chrono::nanoseconds (bug in the implementation of the norm).
    typedef std::chrono::nanoseconds nanoseconds;
    nanoseconds ns (nanoseconds::rep (1000000000.0 * seconds));
    std::this_thread::sleep_for(ns);
  }
} // cpp11

} //namespace CGAL

#endif // CGAL_THREAD_H
