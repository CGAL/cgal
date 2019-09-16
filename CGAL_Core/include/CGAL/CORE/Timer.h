/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 * You can redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: Timer.h
 * Synopsis:
 *      Timer is a class to provide simple timing functions:
 *
 *      Here is an example of how to use it:
 *
 *		Timer timer;
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
 * SPDX-License-Identifier: LGPL-3.0+
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
