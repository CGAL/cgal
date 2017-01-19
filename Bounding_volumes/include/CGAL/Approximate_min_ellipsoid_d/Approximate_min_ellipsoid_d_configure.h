// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

#ifndef CGAL_APPROX_MIN_ELL_CONFIGURE_H
#define CGAL_APPROX_MIN_ELL_CONFIGURE_H

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/basic.h>
#include <iostream>
#include <iomanip>

#if (defined(CGAL_NO_ASSERTIONS) || defined(NDEBUG))
  #undef CGAL_APPEL_ASSERTION_MODE
  #undef CGAL_APPEL_EXP_ASSERTION_MODE
  #undef CGAL_APPEL_LOG_MODE
  #undef CGAL_APPEL_TIMER_MODE
  #undef CGAL_APPEL_STATS_MODE
#else
  #if defined(CGAL_CHECK_EXPENSIVE)
    #define CGAL_APPEL_EXP_ASSERTION_MODE
  #endif
  #define CGAL_APPEL_ASSERTION_MODE
  //#define CGAL_APPEL_LOG_MODE // you might want to temporarily uncomment
                                // this for debugging
#endif

// The following "flags" are best defined at compilation time (e.g.,
// by passing -DCGAL_APPEL_LOG_MODE to the (GCC) compiler):
// #define CGAL_APPEL_ASSERTION_MODE     // enables cheap assertions
// #define CGAL_APPEL_EXP_ASSERTION_MODE // enables cheap & expensive assert's
// #define CGAL_APPEL_LOG_MODE           // enables logging
// #define CGAL_APPEL_TIMER_MODE         // enables runtime measurements
// #define CGAL_APPEL_STATS_MODE         // enables statistics

// The following debugging routines are available with this package:
// (Todo: these routines turned out to not be the most useful for
// debugging; one could improve this, but it's only debugging, so who
// cares ...)
//
// - CGAL_APPEL_ASSERT(expr): Asserts that the expression expr evaluates
//   to true.  This only happens when the assertion mode is enabled.  This
//   macro should only be used for cheap assertions (which do not heavily
//   affect the running time).
//
// - CGAL_APPEL_ASSERT_EXPENSIVE(expr): Asserts that the Boolean expression
//   expr evaluates to true.  This only happens when expensive assertion
//   mode is on, and not otherwise.  (The assertion mode alone doesn't
//   enable such assertions.)
//
// - CGAL_APPEL_LOG(channel,expr): Writes to file channel.log the string
//   produced by the expression 's << expr' (where s is the stream for
//   file channel.log).  This only happens when the logging mode is
//   enabled; nothing is done otherwise.
//   Example: CGAL_APPEL_LOG("debug","i is " << i << std::endl)
//
// - CGAL_APPEL_IF_ASSERTIONS(expr): evaluates the expression expr when
//   assertion mode is on; it does nothing otherwise.
//
// - CGAL_APPEL_TIMER_START(timer): Restarts the timer with name timer.
//   Nothing happens when the timer mode is disabled (see
//   CGAL_APPEL_TIMER_MODE above).
//
// - CGAL_APPEL_TIMER_PRINT(channel,timer,msg): Prints the string msg followed
//   by the time elapsed since the last CGAL_APPEL_TIMER_START(timer).  The
//   destination is either channel 'channel' in case logging mode is on,
//   or std::cout otherwise.
//
// - CGAL_APPEL_TIMER_STRING(timer): Returns a string with the number of
//   seconds elapsed since the last CGAL_APPEL_TIMER_START(timer).
//
// - CGAL_APPEL_IF_TIMING(expr): evaluates the expression when timing mode is
//   enabled (see CGAL_APPEL_TIMER_MODE above), or does nothing otherwise.
//
// - CGAL_APPEL_IF_STATS(expr): evaluates the expression when stats mode is
//   enabled (see CGAL_APPEL_STATS_MODE above), or does nothing otherwise.
//

#ifdef CGAL_APPEL_ASSERTION_MODE
  #define CGAL_APPEL_ASSERT(expr) CGAL_assertion(expr)
#else
  #define CGAL_APPEL_ASSERT(expr)
#endif

#ifdef CGAL_APPEL_EXP_ASSERTION_MODE
  #define CGAL_APPEL_ASSERT_EXPENSIVE(expr) CGAL_assertion(expr)
#else
  #define CGAL_APPEL_ASSERT_EXPENSIVE(expr)
#endif

#ifdef CGAL_APPEL_LOG_MODE
  #define CGAL_APPEL_LOG(channel,expr) \
    { \
      std::stringstream s; \
      s << expr; \
      CGAL::Approximate_min_ellipsoid_d_impl:: \
      Logger::instance().log(channel,s.str()); \
    }
#else
  #define CGAL_APPEL_LOG(channel,expr)
#endif

#if defined(CGAL_APPEL_ASSERTION_MODE)||defined(CGAL_APPEL_EXP_ASSERTION_MODE)
  #define CGAL_APPEL_IF_ASSERTIONS(expr) expr
#else
  #define CGAL_APPEL_IF_ASSERTIONS(expr)
#endif

#ifdef CGAL_APPEL_TIMER_MODE
  #define CGAL_APPEL_TIME(expr) expr
  #define CGAL_APPEL_TIMER_START(timer) \
    { \
      CGAL::Approximate_min_ellipsoid_d_impl::Timer::instance().start(timer); \
    }
  #ifdef CGAL_APPEL_LOG_MODE
    #define CGAL_APPEL_TIMER_PRINT(channel,timer,msg) \
      { \
        CGAL_APPEL_LOG(channel,msg \
	  	          << std::setiosflags(std::ios::fixed) \
                          << std::setprecision(5) \
                          << CGAL::Approximate_min_ellipsoid_d_impl:: \
                             Timer::instance().lapse(timer) \
                          << 's' << std::endl) \
      }
   #else
     #define CGAL_APPEL_TIMER_PRINT(channel,timer,msg) \
       { \
         std::cerr << msg << std::setiosflags(std::ios::fixed) \
                   << std::setprecision(5) \
                   << CGAL::Approximate_min_ellipsoid_d_impl:: \
                      Timer::instance().lapse(timer) \
                   << 's' << std::endl; \
       }
    #endif
  #define CGAL_APPEL_TIMER_STRING(timer) \
    CGAL::Approximate_min_ellipsoid_d_impl::Timer::instance().lapse(timer)
#else
  #define CGAL_APPEL_TIME(expr)
  #define CGAL_APPEL_TIMER_START(timer)
  #define CGAL_APPEL_TIMER_PRINT(channel,timer,msg)
  #define CGAL_APPEL_TIMER_STRING(timer)
#endif

#ifdef CGAL_APPEL_STATS_MODE
  #define CGAL_APPEL_IF_STATS(expr) expr
#else
  #define CGAL_APPEL_IF_STATS(expr)
#endif

#include <CGAL/Approximate_min_ellipsoid_d/Approximate_min_ellipsoid_d_debug.h>

#endif // CGAL_APPROX_MIN_ELL_CONFIGURE_H
