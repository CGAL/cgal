//  (C) Copyright John Maddock 2001 - 2003. 
//  (C) Copyright Darin Adler 2001 - 2002. 
//  (C) Copyright Bill Kempf 2002. 
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for most recent version.

//  Mac OS specific config options:

#define BOOST_PLATFORM "Mac OS"

#if __MACH__ && !defined(_MSL_USING_MSL_C)

// Using the Mac OS X system BSD-style C library.

#  ifndef BOOST_HAS_UNISTD_H
#    define BOOST_HAS_UNISTD_H
#  endif
// boilerplate code:
#  ifndef TARGET_CARBON
#     include <boost/config/posix_features.hpp>
#  endif
#  ifndef BOOST_HAS_STDINT_H
#     define BOOST_HAS_STDINT_H
#  endif

//
// BSD runtime has pthreads, sigaction, sched_yield and gettimeofday,
// of these only pthreads are advertised in <unistd.h>, so set the 
// other options explicitly:
//
#  define BOOST_HAS_SCHED_YIELD
#  define BOOST_HAS_GETTIMEOFDAY
#  define BOOST_HAS_SIGACTION

#  if (__GNUC__ < 3) && !defined( __APPLE_CC__)

// GCC strange "ignore std" mode works better if you pretend everything
// is in the std namespace, for the most part.

#    define BOOST_NO_STDC_NAMESPACE
#  endif

#else

// Using the MSL C library.

// We will eventually support threads in non-Carbon builds, but we do
// not support this yet.
#  if TARGET_CARBON

#    define BOOST_HAS_MPTASKS

// The MP task implementation of Boost Threads aims to replace MP-unsafe
// parts of the MSL, so we turn on threads unconditionally.
#    define BOOST_HAS_THREADS

// The remote call manager depends on this.
#    define BOOST_BIND_ENABLE_PASCAL

#  endif

#endif



