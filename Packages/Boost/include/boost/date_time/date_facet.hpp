#ifndef _DATE_TIME_DATE_FACET__HPP___
#define _DATE_TIME_DATE_FACET__HPP___

/* Copyright (c) 2004 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author:  Martin Andrian, Jeff Garland
 * $Date$
 */


#include "boost/algorithm/string.hpp" //todo narrow this
#include "boost/date_time/period.hpp"
#include "boost/date_time/special_values_formatter.hpp"
#include "boost/date_time/period_formatter.hpp"
#include <string>
#include <vector>

namespace boost { namespace date_time {


  template <class date_type,
            class CharT, 
            class OutItrT = std::ostreambuf_iterator<CharT, std::char_traits<CharT> > >
  class date_facet : public std::locale::facet {
  public:
    typedef typename date_type::duration_type duration_type;
    typedef boost::date_time::period<date_type,duration_type> period_type;
    typedef std::basic_string<CharT> string_type;
    typedef CharT                    char_type;
    typedef boost::date_time::period_formatter<CharT>  period_formatter_type;
    typedef boost::date_time::special_values_formatter<CharT>  special_values_formatter_type;
    typedef std::vector<std::basic_string<CharT> > input_collection_type;
    static const char_type long_weekday_format[3];
    static const char_type short_weekday_format[3];
    static const char_type long_month_format[3];
    static const char_type short_month_format[3];
    static const char_type default_period_separator[4];
    static const char_type standard_format_specifier[3];
    static const char_type iso_format_specifier[7];
    static const char_type iso_format_extended_specifier[9];
    static std::locale::id id;
    
    explicit date_facet(::size_t a_ref = 0) 
      : std::locale::facet(a_ref), 
        m_format(standard_format_specifier)
    {}

    explicit date_facet(const char_type* format,
                        const input_collection_type& short_month_names,
                        ::size_t ref_count = 0) 
      : std::locale::facet(ref_count), 
        m_format(format),
        m_month_short_names(short_month_names) 
    {}

     
    explicit date_facet(const char_type* format, 
                        period_formatter_type period_formatter = period_formatter_type(), 
                        special_values_formatter_type special_values_formatter = special_values_formatter_type(),
                        ::size_t ref_count = 0)
      : std::locale::facet(ref_count), 
        m_format(format), 
        m_period_formatter(period_formatter),
        m_special_values_formatter(special_values_formatter)
     {}
    void format(const char_type* const format) {
      m_format = format;
    }
    virtual void set_iso_format()
    {
      m_format = iso_format_specifier;
    }
    virtual void set_iso_extended_format()
    {
      m_format = iso_format_extended_specifier;
    }
    
    void period_formatter(period_formatter_type period_formatter) {
      m_period_formatter= period_formatter;
    }
    void special_value_formatting(const special_values_formatter_type& svf) 
    {
      m_special_values_formatter = svf;
    }
    void short_weekday_names(const input_collection_type& short_weekday_names)
    {
      m_weekday_short_names = short_weekday_names;
    }
    void long_weekday_names(const input_collection_type& long_weekday_names)
    {
      m_weekday_long_names = long_weekday_names;
    }

    void short_month_names(const input_collection_type& short_month_names)
    {
      m_month_short_names = short_month_names;
    }

    void long_month_names(const input_collection_type& long_month_names)
    {
      m_month_long_names = long_month_names;
    }

    OutItrT put(OutItrT next, 
                std::ios_base& a_ios, 
                char_type fill_char, 
                const date_type& d) const 
    {
      if (d.is_special()) { 
        return do_put_special(next, a_ios, fill_char, d.as_special());
      }
      //The following line of code required the date to support a to_tm function
      return do_put_tm(next, a_ios, fill_char, to_tm(d), m_format);
    }

    // date durations don't interact with the facet...no put code needed here.
    // todo is that really true?  What about special values?

    OutItrT put(OutItrT next, 
                std::ios_base& a_ios, 
                char_type fill_char, 
                const period_type& p) const 
    {
      return m_period_formatter.put_period(next, a_ios, fill_char, p, *this);
    }
  protected:
    virtual OutItrT do_put_special(OutItrT next, 
                                   std::ios_base& a_ios, 
                                   char_type fill_char, 
                                   const boost::date_time::special_values sv) const 
    {
      m_special_values_formatter.put_special(next, sv);
      return next;
    }
    virtual OutItrT do_put_tm(OutItrT next, 
                              std::ios_base& a_ios, 
                              char_type fill_char, 
                              const tm& tm_value,
                              string_type a_format) const 
    {
      // update format string with custom names
      if (m_weekday_long_names.size()) {
        boost::algorithm::replace_all(a_format, 
                                      long_weekday_format, 
                                      m_weekday_long_names[tm_value.tm_wday]);
      }
      if (m_weekday_short_names.size()) {
        boost::algorithm::replace_all(a_format, 
                                      short_weekday_format, 
                                      m_weekday_short_names[tm_value.tm_wday]);

      }
      if (m_month_long_names.size()) {
        boost::algorithm::replace_all(a_format, 
                                      long_month_format, 
                                      m_month_long_names[tm_value.tm_mon]);
      }
      if (m_month_short_names.size()) {
        boost::algorithm::replace_all(a_format, 
                                      short_month_format, 
                                      m_month_short_names[tm_value.tm_mon]);
      }
      // use time_put facet to create final string
      return std::use_facet<std::time_put<CharT> >(a_ios.getloc()).put(next, a_ios, 
                                                                       fill_char, 
                                                                       &tm_value,
                                                                       &*a_format.begin(), 
                                                                       &*a_format.end());
    }
  protected:
    string_type                   m_format;
    period_formatter_type         m_period_formatter;
    special_values_formatter_type m_special_values_formatter;
    input_collection_type         m_month_short_names;
    input_collection_type         m_month_long_names;
    input_collection_type         m_weekday_short_names;
    input_collection_type         m_weekday_long_names;
  private:
  };

  template <class date_type, class CharT, class OutItrT>
  std::locale::id date_facet<date_type, CharT, OutItrT>::id;

  template <class date_type, class CharT, class OutItrT>  
  const typename date_facet<date_type, CharT, OutItrT>::char_type 
  date_facet<date_type, CharT, OutItrT>::long_weekday_format[3] = {'%','A'};

  template <class date_type, class CharT, class OutItrT>  
  const typename date_facet<date_type, CharT, OutItrT>::char_type 
  date_facet<date_type, CharT, OutItrT>::short_weekday_format[3] = {'%','a'};

  template <class date_type, class CharT, class OutItrT>  
  const typename date_facet<date_type, CharT, OutItrT>::char_type 
  date_facet<date_type, CharT, OutItrT>::long_month_format[3] = {'%','B'};

  template <class date_type, class CharT, class OutItrT>  
  const typename date_facet<date_type, CharT, OutItrT>::char_type 
  date_facet<date_type, CharT, OutItrT>::short_month_format[3] = {'%','b'};

  template <class date_type, class CharT, class OutItrT>  
  const typename date_facet<date_type, CharT, OutItrT>::char_type 
  date_facet<date_type, CharT, OutItrT>::default_period_separator[4] = { ' ', '/', ' '};

  template <class date_type, class CharT, class OutItrT>  
  const typename date_facet<date_type, CharT, OutItrT>::char_type 
  date_facet<date_type, CharT, OutItrT>::standard_format_specifier[3] = 
    {'%', 'x' };

  template <class date_type, class CharT, class OutItrT>  
  const typename date_facet<date_type, CharT, OutItrT>::char_type 
  date_facet<date_type, CharT, OutItrT>::iso_format_specifier[7] = 
    {'%', 'Y', '%', 'm', '%', 'd' };

  template <class date_type, class CharT, class OutItrT>  
  const typename date_facet<date_type, CharT, OutItrT>::char_type 
  date_facet<date_type, CharT, OutItrT>::iso_format_extended_specifier[9] = 
    {'%', 'Y', '-', '%', 'm', '-', '%', 'd' };

} }


#endif
