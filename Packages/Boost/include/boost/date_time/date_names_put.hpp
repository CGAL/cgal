#ifndef DATE_TIME_DATE_NAMES_PUT_HPP___
#define DATE_TIME_DATE_NAMES_PUT_HPP___

/* Copyright (c) 2002,2003 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland 
 * $Date$
 */


#include "boost/date_time/locale_config.hpp" // set BOOST_DATE_TIME_NO_LOCALE

#ifndef BOOST_DATE_TIME_NO_LOCALE

#include "boost/date_time/special_defs.hpp"
#include "boost/date_time/date_defs.hpp"
#include "boost/date_time/parse_format_base.hpp"
#include <locale>


namespace boost {
namespace date_time {

    //! Output facet base class for gregorian dates.
    /*! This class is a base class for date facets used to localize the
     *  names of months and the names of days in the week.
     * 
     * Requirements of Config
     *  - define an enumeration month_enum that enumerates the months. 
     *    The enumeration should be '1' based eg: Jan==1
     *  - define as_short_string and as_long_string
     *
     * (see langer & kreft p334).
     * 
     */
    template<class Config,
             class charT = char, 
             class OutputIterator = std::ostreambuf_iterator<charT> >
    class date_names_put : public std::locale::facet
    {
    public:
      date_names_put() {};
      typedef OutputIterator iter_type;
      typedef typename Config::month_type month_type;
      typedef typename Config::month_enum month_enum;
      typedef typename Config::weekday_enum weekday_enum;
      typedef typename Config::special_value_enum special_value_enum;
      //typedef typename Config::format_type format_type;
      typedef std::basic_string<charT> string_type;

      static std::locale::id id;
      void put_special_value(iter_type& oitr, special_value_enum sv) const
      {
        do_put_special_value(oitr, sv);
      }
      void put_month_short(iter_type& oitr, month_enum moy) const
      {
        do_put_month_short(oitr, moy);
      }
      void put_month_long(iter_type& oitr, month_enum moy) const
      {
        do_put_month_long(oitr, moy);
      }
      void put_weekday_short(iter_type& oitr, weekday_enum wd) const
      {
        do_put_weekday_short(oitr, wd);
      }
      void put_weekday_long(iter_type& oitr, weekday_enum wd) const
      {
        do_put_weekday_long(oitr, wd);
      }
      bool has_date_sep_chars() const
      {
        return do_has_date_sep_chars();
      }
      void year_sep_char(iter_type& oitr) const
      {
        do_year_sep_char(oitr);
      }
      //! char between year-month
      void month_sep_char(iter_type& oitr) const
      {
        do_month_sep_char(oitr);
      }
      //! Char to separate month-day
      void day_sep_char(iter_type& oitr) const
      {
        do_day_sep_char(oitr);
      }
      //! Determines the order to put the date elements
      ymd_order_spec date_order() const
      {
        return do_date_order();
      }
      //! Determines if month is displayed as integer, short or long string
      month_format_spec month_format() const
      {
        return do_month_format();
      }

    protected:
      //! Default facet implementation uses month_type defaults
      virtual void do_put_month_short(iter_type& oitr, month_enum moy) const
      {
        month_type gm(moy);
        put_string(oitr, gm.as_short_string());
      }
      //! Default facet implementation uses month_type defaults
      virtual void do_put_month_long(iter_type& oitr, 
                                     month_enum moy) const
      {
        month_type gm(moy);
        put_string(oitr, gm.as_long_string());
      }
      //! Default facet implementation for special value types
      virtual void do_put_special_value(iter_type& oitr, special_value_enum sv) const
      {
        switch (sv) {
          case not_a_date_time: 
          { 
            put_string(oitr, "not-a-date-time");
            break;
          }
          case pos_infin: 
          { 
            put_string(oitr, "+infinity");
            break;
          }
          case neg_infin: 
          { 
            put_string(oitr, "-infinity");
            break;
          }
          default: {} //quiet compilers that want all cases covered here (eg: gcc 3.1)
        }
      }
      virtual void do_put_weekday_short(iter_type& oitr, weekday_enum wd) const
      {
      }
      virtual void do_put_weekday_long(iter_type& oitr, weekday_enum wd) const
      {
      }
      virtual bool do_has_date_sep_chars() const
      {
        return true;
      }
      virtual void do_year_sep_char(iter_type& oitr) const
      {
        put_string(oitr, "-");
      }
      //! char between year-month
      virtual void do_month_sep_char(iter_type& oitr) const
      {
        put_string(oitr, "-");
      }
      //! Char to separate month-day
      virtual void do_day_sep_char(iter_type& oitr) const
      {
        put_string(oitr, "-");
      }
      //! Default for date order 
      virtual ymd_order_spec do_date_order() const
      {
        return ymd_order_iso;
      }
      //! Default month format
      virtual month_format_spec do_month_format() const
      {
        return month_as_short_string;
      }
      void put_string(iter_type& oi, const char* const s) const
      {
        string_type s1(s);
        typename string_type::iterator si,end;
        for (si=s1.begin(), end=s1.end(); si!=end; si++, oi++) {
          *oi = *si;
        }
      }
    };
    
    //! Generate storage location for a std::locale::id 
    template<class Config, class charT, class OutputIterator>
    std::locale::id date_names_put<Config, charT, OutputIterator>::id;

    //! An date name output facet that takes an array of char* to define strings
    template<class Config,
             class charT = char, 
             class OutputIterator = std::ostreambuf_iterator<charT> >
    class all_date_names_put : public date_names_put<Config, charT, OutputIterator>
    {
    public:
      all_date_names_put(const char* const month_short_names[],
                         const char* const month_long_names[],
                         const char* const special_value_names[],
                         const char* const weekday_short_names[],
                         const char* const weekday_long_names[],
                         char separator_char = '-',
                         ymd_order_spec order_spec = ymd_order_iso,
                         month_format_spec month_format = month_as_short_string) :
        month_short_names_(month_short_names),
        month_long_names_(month_long_names),
        special_value_names_(special_value_names),
        weekday_short_names_(weekday_short_names),
        weekday_long_names_(weekday_long_names),
        order_spec_(order_spec),
        month_format_spec_(month_format)
      {
        separator_char_[0] = separator_char;
        separator_char_[1] = '\0';

      };
      typedef OutputIterator iter_type;
      typedef typename Config::month_enum month_enum;
      typedef typename Config::weekday_enum weekday_enum;
      typedef typename Config::special_value_enum special_value_enum;

    protected:
      //! Generic facet that takes array of chars
      virtual void do_put_month_short(iter_type& oitr, month_enum moy) const
      {
        this->put_string(oitr, month_short_names_[moy-1]);
      }
      //! Long month names 
      virtual void do_put_month_long(iter_type& oitr, month_enum moy) const
      {
        this->put_string(oitr, month_long_names_[moy-1]);
      }
      //! Special values names
      virtual void do_put_special_value(iter_type& oitr, special_value_enum sv) const
      {
        this->put_string(oitr, special_value_names_[sv]);
      }
      virtual void do_put_weekday_short(iter_type& oitr, weekday_enum wd) const
      {
        this->put_string(oitr, weekday_short_names_[wd]);
      }
      virtual void do_put_weekday_long(iter_type& oitr, weekday_enum wd) const
      {
        this->put_string(oitr, weekday_long_names_[wd]);
      }
      //! char between year-month
      virtual void do_month_sep_char(iter_type& oitr) const
      {
        this->put_string(oitr, separator_char_);
      }
      //! Char to separate month-day
      virtual void do_day_sep_char(iter_type& oitr) const
      {
        this->put_string(oitr, separator_char_);
      }
      //! Set the date ordering
      virtual ymd_order_spec do_date_order() const
      {
        return order_spec_;
      }
      //! Set the date ordering
      virtual month_format_spec do_month_format() const
      {
        return month_format_spec_;
      }

    private:
      const char* const* month_short_names_;
      const char* const* month_long_names_;
      const char* const* special_value_names_;
      const char* const* weekday_short_names_;
      const char* const* weekday_long_names_;
      char separator_char_[2];
      ymd_order_spec order_spec_;
      month_format_spec month_format_spec_;      
    };

} } //namespace boost::date_time

#endif //BOOST_NO_STD_LOCALE

#endif
