#ifndef DATE_TIME_C_TIME_HPP___
#define DATE_TIME_C_TIME_HPP___

/* Copyright (c) 2002,2003 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland 
 * $Date$
 */


/*! @file c_time.hpp
  Provide workarounds related to the ctime header
*/

#include "boost/date_time/compiler_config.hpp"
#include <ctime>
//Work around libraries that don't put time_t and time in namespace std

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std { using ::time_t; using ::time; using ::localtime;
                using ::tm;  using ::gmtime; }
#endif

//The following is used to support high precision time clocks
#ifdef BOOST_HAS_GETTIMEOFDAY
#include <sys/time.h>
#endif


#endif
