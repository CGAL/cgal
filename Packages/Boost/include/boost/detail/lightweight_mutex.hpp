#ifndef BOOST_DETAIL_LIGHTWEIGHT_MUTEX_HPP_INCLUDED
#define BOOST_DETAIL_LIGHTWEIGHT_MUTEX_HPP_INCLUDED

// MS compatible compilers support #pragma once

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

//
//  boost/detail/lightweight_mutex.hpp - lightweight mutex
//
//  Copyright (c) 2002, 2003 Peter Dimov and Multi Media Ltd.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//  typedef <unspecified> boost::detail::lightweight_mutex;
//
//  boost::detail::lightweight_mutex meets a subset of the Mutex concept
//  requirements: http://www.boost.org/libs/thread/doc/mutex_concept.html#Mutex
//
//  * Used by the smart pointer library
//  * Performance oriented
//  * Header-only implementation
//  * Small memory footprint
//  * Not a general purpose mutex, use boost::mutex, CRITICAL_SECTION or
//    pthread_mutex instead.
//  * Never spin in a tight lock/do-something/unlock loop, since
//    lightweight_mutex does not guarantee fairness.
//  * Never keep a lightweight_mutex locked for long periods.
//
//  The current implementation can use a pthread_mutex, a CRITICAL_SECTION,
//  or a platform-specific spinlock.
//
//  You can force a particular implementation by defining BOOST_LWM_USE_PTHREADS,
//  BOOST_LWM_USE_CRITICAL_SECTION, or BOOST_LWM_USE_SPINLOCK.
//
//  If neither macro has been defined, the default is to use a spinlock on Win32,
//  and a pthread_mutex otherwise.
//
//  Note that a spinlock is not a general synchronization primitive. In particular,
//  it is not guaranteed to be a memory barrier, and it is possible to "livelock"
//  if a lower-priority thread has acquired the spinlock but a higher-priority
//  thread is spinning trying to acquire the same lock.
//
//  For these reasons, spinlocks have been disabled by default except on Windows,
//  where a spinlock can be several orders of magnitude faster than a CRITICAL_SECTION.


//  Note: lwm_linux.hpp has been disabled by default; see the comments
//        inside for more info.


#include <boost/config.hpp>

//  Note to implementors: if you write a platform-specific spinlock
//  for a platform that supports pthreads, be sure to test its performance
//  against the pthreads-based version using shared_ptr_timing_test.cpp and
//  shared_ptr_mt_test.cpp. Custom versions are usually not worth the trouble
//  _unless_ the performance gains are substantial.
//
//  Be sure to compare against a "real" pthreads library;
//  shared_ptr_timing_test.cpp will compile succesfully with a stub do-nothing
//  pthreads library, since it doesn't create any threads.

#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32__)) && !defined(BOOST_LWM_USE_CRITICAL_SECTION) && !defined(BOOST_LWM_USE_PTHREADS)
# define BOOST_LWM_WIN32
#endif

#if !defined(BOOST_HAS_THREADS)
#  if defined(BOOST_LWM_WIN32)
#    include <boost/detail/lwm_win32_nt.hpp>
#  else
#    include <boost/detail/lwm_nop.hpp>
#  endif
#elif defined(BOOST_LWM_USE_SPINLOCK) && defined(BOOST_USE_ASM_ATOMIC_H)
#  include <boost/detail/lwm_linux.hpp>
#elif defined(BOOST_LWM_USE_CRITICAL_SECTION)
#  include <boost/detail/lwm_win32_cs.hpp>
#elif defined(BOOST_LWM_USE_PTHREADS)
#  include <boost/detail/lwm_pthreads.hpp>
#elif defined(BOOST_LWM_WIN32)
#  include <boost/detail/lwm_win32.hpp>
#elif defined(BOOST_LWM_USE_SPINLOCK) && defined(__sgi)
#  include <boost/detail/lwm_irix.hpp>
#elif defined(BOOST_LWM_USE_SPINLOCK) && defined(__GLIBCPP__)
#  include <boost/detail/lwm_gcc.hpp>
#elif defined(BOOST_HAS_PTHREADS)
#  define BOOST_LWM_USE_PTHREADS
#  include <boost/detail/lwm_pthreads.hpp>
#else
// Use #define BOOST_DISABLE_THREADS to avoid the error
#  error Unrecognized threading platform
#endif

#endif // #ifndef BOOST_DETAIL_LIGHTWEIGHT_MUTEX_HPP_INCLUDED
