#ifndef POSIX_TIME_CONVERSION_HPP___
#define POSIX_TIME_CONVERSION_HPP___

/* Copyright (c) 2002,2003 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland 
 * $Date$
 */

#include "boost/date_time/posix_time/ptime.hpp"
#include "boost/date_time/posix_time/posix_time_duration.hpp"
#include "boost/date_time/c_time.hpp"

namespace boost {

namespace posix_time {


  //! Function that converts a time_t into a ptime.
  inline
  ptime from_time_t(std::time_t t) 
  {
    ptime start(gregorian::date(1970,1,1));
    return start + seconds(t);
  }


} } //namespace boost::posix_time




#endif

