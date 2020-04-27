/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: Timer.h
 * Synopsis:
 *      Timer is a class to provide simple timing functions:
 *
 *      Here is an example of how to use it:
 *
 *                Timer timer;
 *
 *              timer.start();
 *              .. do some tasks for timing ..
 *              timer.stop();
 *
 *              long clock = timer.getClocks();     // get CPU clocks
 *              long seconds = time.getSeconds();   // get seconds
 *
 * Written by
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/

#ifndef _CORE_TIMER_H_
#define _CORE_TIMER_H_

#include <CGAL/CORE/Impl.h>
#include <ctime>

namespace CORE {

class Timer {
private:
  long startClock;
  long clocks;

public:
  Timer() : startClock(0), clocks(0) {}

  void start() {
    startClock = clock();
  }

  void stop() {
    clocks = clock() - startClock;
  }

  long getClocks() {
    return clocks;
  }

  float getSeconds() {
    return (float)clocks / CLOCKS_PER_SEC;
  }
};

} //namespace CORE
#endif // _CORE_TIMER_H_
