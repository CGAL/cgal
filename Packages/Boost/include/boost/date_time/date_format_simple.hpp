#ifndef DATE_TIME_SIMPLE_FORMAT_HPP___
#define DATE_TIME_SIMPLE_FORMAT_HPP___

/* Copyright (c) 2002,2003 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst
 * $Date$
 */

#include "boost/date_time/parse_format_base.hpp"

namespace boost {
namespace date_time {

//! Class to provide simple basic formatting rules
class simple_format {
public:

  //! String used printed is date is invalid
  static const char* not_a_date()
  {
    return "not-a-date-time";
  }
  //! String used to for positive infinity value
  static const char* pos_infinity()
  {       //2001-Jan-03
    return "+infinity";
  }
  //! String used to for positive infinity value
  static const char* neg_infinity()
  {
    return "-infinity";
  }
  //! Describe month format
  static month_format_spec month_format()
  {
    return month_as_short_string;
  }
  static ymd_order_spec date_order()
  {
    return ymd_order_iso; //YYYY-MM-DD
  }
  //! This format uses '-' to separate date elements
  static bool has_date_sep_chars()
  {
    return true;
  }
  //! Char to sep?
  static char year_sep_char()
  {
    return '-';
  }
  //! char between year-month
  static char month_sep_char()
  {
    return '-';
  }
  //! Char to separate month-day
  static char day_sep_char()
  {
    return '-';
  }
  //! char between date-hours
  static char hour_sep_char()
  {
    return ' ';
  }
  //! char between hour and minute
  static char minute_sep_char()
  {
    return ':';
  }
  //! char for second
  static char second_sep_char()
  {
    return ':';
  }

};

} } //namespace date_time




#endif
