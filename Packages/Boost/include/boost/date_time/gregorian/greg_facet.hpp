#ifndef GREGORIAN_FACET_HPP___
#define GREGORIAN_FACET_HPP___

/* Copyright (c) 2002,2003 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst
 * $Date$
 */

#include "boost/date_time/gregorian/gregorian_types.hpp"
#include "boost/date_time/date_formatting_locales.hpp" // sets BOOST_DATE_TIME_NO_LOCALE
#include <string>

//This file is basically commented out if locales are not supported
#ifndef BOOST_DATE_TIME_NO_LOCALE


namespace boost {
namespace gregorian {

  //! Configuration of the output facet template
  struct greg_facet_config
  {
    typedef boost::gregorian::greg_month month_type;
    typedef boost::date_time::special_values special_value_enum;
    typedef boost::gregorian::months_of_year month_enum;
    typedef boost::date_time::weekdays weekday_enum;
  };

  //! Create the base facet type for gregorian::date
  typedef boost::date_time::date_names_put<greg_facet_config> greg_base_facet;
  
  //! ostream operator for gregorian::date
  /*! Uses the date facet to determine various output parameters including:
   *  - string values for the month (eg: Jan, Feb, Mar) (default: English)
   *  - string values for special values (eg: not-a-date-time) (default: English)
   *  - selection of long, short strings, or numerical month representation (default: short string)
   *  - month day year order (default yyyy-mmm-dd)
   */
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, const date& d)
  {
    typedef boost::date_time::ostream_date_formatter<date, greg_base_facet, charT> greg_ostream_formatter;
    greg_ostream_formatter::date_put(d, os);
    return os;
  }

  //! operator<< for gregorian::greg_month typically streaming: Jan, Feb, Mar...
  /*! Uses the date facet to determine output string as well as selection of long or short strings.
   *  Default if no facet is installed is to output a 2 wide numeric value for the month
   *  eg: 01 == Jan, 02 == Feb, ... 12 == Dec.
   */
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, const greg_month& m)
  {
    typedef boost::date_time::ostream_month_formatter<greg_base_facet, charT> greg_month_formatter;
    std::locale locale = os.getloc();
    if (std::has_facet<greg_base_facet>(locale)) {
      const greg_base_facet& f = std::use_facet<greg_base_facet>(locale);
      greg_month_formatter::format_month(m, os, f);

    }
    else { //default to numeric
      os  << std::setw(2) << std::setfill('0') << m.as_number();
    }

    return os;
  }

  //! operator<< for gregorian::greg_weekday typically streaming: Sun, Mon, Tue, ...
  /*! Uses the date facet to determine output string as well as selection of long or short string.
   *  Default if no facet is installed is to output a 3 char english string for the
   *  day of the week.
   */
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, const greg_weekday& wd)
  {
    typedef boost::date_time::ostream_weekday_formatter<greg_weekday, greg_base_facet, charT> greg_weekday_formatter;
    std::locale locale = os.getloc();
    if (std::has_facet<greg_base_facet>(locale)) {
      const greg_base_facet& f = std::use_facet<greg_base_facet>(locale);
      greg_weekday_formatter::format_weekday(wd.as_enum(), os, f, true);
    }
    else { //default to short English string eg: Sun, Mon, Tue, Wed...
      os  << wd.as_short_string();
    }

    return os;
  }

  //! operator<< for gregorian::date_period typical output: [2002-Jan-01/2002-Jan-31]
  /*! Uses the date facet to determine output string as well as selection of long 
   *  or short string fr dates.
   *  Default if no facet is installed is to output a 3 char english string for the
   *  day of the week.
   */
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, const date_period& dp)
  {
    os << '['; //TODO: facet or manipulator for periods?
    os << dp.begin();
    os << '/'; //TODO: facet or manipulator for periods?
    os << dp.last();
    os << ']'; 
    return os;
  }

  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, const date_duration& dd)
  {
    //os << dd.days();
    os << dd.get_rep();
    return os;
  }

  //! operator<< for gregorian::partial_date. Output: "Jan 1"
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, const partial_date& pd)
  {
    os << pd.day() << ' ' << pd.month().as_short_string() ; 
    return os;
  }

  //! operator<< for gregorian::nth_kday_of_month. Output: "first Mon of Jun"
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, 
             const nth_kday_of_month& nkd)
  {
    os << nkd.nth_week_as_str() << ' ' 
       << nkd.day_of_week() << " of "
       << nkd.month().as_short_string() ; 
    return os;
  }

  //! operator<< for gregorian::first_kday_of_month. Output: "first Mon of Jun"
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, 
             const first_kday_of_month& fkd)
  {
    os << "first " << fkd.day_of_week() << " of " 
       << fkd.month().as_short_string() ; 
    return os;
  }

  //! operator<< for gregorian::last_kday_of_month. Output: "last Mon of Jun"
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, 
             const last_kday_of_month& lkd)
  {
    os << "last " << lkd.day_of_week() << " of " 
       << lkd.month().as_short_string() ; 
    return os;
  }

  //! operator<< for gregorian::first_kday_after. Output: "first Mon after"
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, 
             const first_kday_after& fka)
  {
    os << fka.day_of_week() << " after"; 
    return os;
  }

  //! operator<< for gregorian::first_kday_before. Output: "first Mon before"
  template <class charT, class traits>
  inline
  std::basic_ostream<charT, traits>&
  operator<<(std::basic_ostream<charT, traits>& os, 
             const first_kday_before& fkb)
  {
    os << fkb.day_of_week() << " before"; 
    return os;
  }

} } //namespace gregorian

#endif  
    
    
#endif

