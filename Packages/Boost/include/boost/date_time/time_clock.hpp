#ifndef DATE_TIME_TIME_CLOCK_HPP___
#define DATE_TIME_TIME_CLOCK_HPP___

/* Copyright (c) 2002,2003 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst 
 * $Date$
 */

/*! @file time_clock.hpp
  This file contains the interface for clock devices.  
*/

#include "boost/date_time/c_time.hpp"

namespace boost {
namespace date_time {


  //! A clock providing time level services based on C time_t capabilities
  /*! This clock provides resolution to the 1 second level
   */
  template<class date_type, class time_type> 
  class second_clock
  {
  public:
    //    typedef typename time_type::date_type date_type;
    typedef typename time_type::time_duration_type time_duration_type;

    static time_type local_time() 
    {
      ::std::time_t t;
      ::std::time(&t); 
      ::std::tm* curr = ::std::localtime(&t);
      return create_time(curr);
    }

    //! Get the current day in universal date as a ymd_type
    static time_type universal_time() 
    {

      ::std::time_t t;
      ::std::time(&t);
      ::std::tm* curr= ::std::gmtime(&t);
      return create_time(curr);
    }

  private:
    static time_type create_time(::std::tm* current)
    {
      date_type d(static_cast<unsigned short>(current->tm_year + 1900), 
                  static_cast<unsigned short>(current->tm_mon + 1), 
                  static_cast<unsigned short>(current->tm_mday));
      time_duration_type td(current->tm_hour,
                            current->tm_min,
                            current->tm_sec);
      return time_type(d,td);
    }
    
  };


} } //namespace date_time


#endif
