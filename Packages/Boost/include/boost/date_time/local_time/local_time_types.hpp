#ifndef LOCAL_TIME_LOCAL_TIME_TYPES_HPP__
#define LOCAL_TIME_LOCAL_TIME_TYPES_HPP__

/* Copyright (c) 2003-2004 CrystalClear Software, Inc.
 * Subject to the Boost Software License, Version 1.0. 
 * (See accompanying file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst
 * $Date$
 */

#include "boost/date_time/local_time/local_date_time.hpp"
#include "boost/date_time/period.hpp"
#include "boost/date_time/time_iterator.hpp"
#include "boost/date_time/compiler_config.hpp"
#if defined(BOOST_DATE_TIME_OPTIONAL_GREGORIAN_TYPES)
#include "boost/date_time/local_time/date_duration_operators.hpp"
#endif //BOOST_DATE_TIME_OPTIONAL_GREGORIAN_TYPES
#include "boost/date_time/local_time/custom_time_zone.hpp"

namespace boost {
namespace local_time {

  typedef boost::date_time::period<local_date_time, 
                                   boost::posix_time::time_duration> local_time_period;

  typedef date_time::time_itr<local_date_time> local_time_iterator;

  typedef date_time::second_clock<local_date_time> local_sec_clock; 
  typedef date_time::microsec_clock<local_date_time> local_microsec_clock;
  //! Shared Pointer for time_zone and posix_time_zone objects
  typedef boost::shared_ptr<time_zone> time_zone_ptr;
  
  //bring special enum values into the namespace
  using date_time::special_values;
  using date_time::not_special;
  using date_time::neg_infin;
  using date_time::pos_infin;
  using date_time::not_a_date_time;
  using date_time::max_date_time;
  using date_time::min_date_time;

}} // namespaces

#endif // LOCAL_TIME_LOCAL_TIME_TYPES_HPP__
