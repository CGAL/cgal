
#ifndef POSIX_TIME_DATE_TIME_FACET__HPP__
#define POSIX_TIME_DATE_TIME_FACET__HPP__

/* Copyright (c) 2004 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author:  Martin Andrian, Jeff Garland
 * $Date$
 */

#include "boost/date_time/date_facet.hpp"

namespace boost {
namespace date_time {

  template <class time_type,
            class CharT, 
            class OutItrT = std::ostreambuf_iterator<CharT, std::char_traits<CharT> > >
  class time_facet : 
    public boost::date_time::date_facet<typename time_type::date_type , CharT, OutItrT> {
   public:
    typedef typename time_type::date_type date_type;
    typedef typename time_type::time_duration_type time_duration_type;
    typedef boost::date_time::period<time_type,time_duration_type> period_type;
    typedef boost::date_time::date_facet<boost::gregorian::date, CharT, OutItrT> base_type;
    typedef typename base_type::string_type string_type;
    typedef typename base_type::char_type   char_type;
    typedef typename base_type::period_formatter_type period_formatter_type;
    typedef typename base_type::special_values_formatter_type special_values_formatter_type;
    static const char_type fractional_seconds_format[3];
    static const char_type fractional_seconds_or_none_format[3];
    static const char_type seconds_with_fractional_seconds_format[3];
    static const char_type seconds_format[3];
    static const char_type standard_format[3];
    static const char_type duration_seperator[2];
    static const char_type iso_time_format_specifier[16];
    static const char_type iso_time_format_extended_specifier[20];
    static std::locale::id id;
    

    explicit time_facet(::size_t a_ref = 0) 
      : base_type(standard_format),
        m_time_duration_seperator(duration_seperator)

    {
      set_duration_seperator_from_facet(); //todo fix this
    }

    explicit time_facet(const char_type* a_format,
                        period_formatter_type period_formatter = period_formatter_type(), 
                        const special_values_formatter_type& special_value_formatter = special_values_formatter_type(), 
                         ::size_t a_ref = 0) 
      : base_type(a_format, 
                  period_formatter,
                  special_value_formatter, 
                  a_ref),
        m_time_duration_seperator(duration_seperator)
    {}


    virtual void set_iso_format()
    {
      m_format = iso_time_format_specifier;
    }
    virtual void set_iso_extended_format()
    {
      m_format = iso_time_format_extended_specifier;
    }

    OutItrT put(OutItrT a_next, 
                std::ios_base& a_ios, 
                char_type a_fill, 
                const time_type& a_time) const 
    {
      if (a_time.is_special()) { 
        return do_put_special(a_next, a_ios, a_fill, 
                              a_time.date().as_special());
      }
      string_type format(m_format);
      string_type frac_str;
      if (format.find(seconds_with_fractional_seconds_format)) {
        // replace %s with %S.nnn 
        frac_str = 
          fractional_seconds_as_string(a_time.time_of_day(), false);
        char_type sep = std::use_facet<std::numpunct<char_type> >(a_ios.getloc()).decimal_point();
        
        string_type replace_string(seconds_format);
        replace_string += sep;
        replace_string += frac_str;
        boost::algorithm::replace_all(format, 
                                      seconds_with_fractional_seconds_format, 
                                      replace_string);
      }
      if (format.find(fractional_seconds_format)) {
        // replace %f with nnnnnnn
        if (!frac_str.size()) {
          frac_str = fractional_seconds_as_string(a_time.time_of_day(), false);
        }
        boost::algorithm::replace_all(format,
                                      fractional_seconds_format, 
                                      frac_str);
      }

      if (format.find(fractional_seconds_or_none_format)) {
        // replace %F with nnnnnnn or nothing if fs == 0
        frac_str = 
          fractional_seconds_as_string(a_time.time_of_day(), true);
        if (frac_str.size()) {
          char_type sep = std::use_facet<std::numpunct<char_type> >(a_ios.getloc()).decimal_point();
          string_type replace_string;
          replace_string += sep;
          replace_string += frac_str;
          boost::algorithm::replace_all(format,
                                        fractional_seconds_or_none_format, 
                                        replace_string);
        }
        else {
          boost::algorithm::erase_all(format,
                                      fractional_seconds_or_none_format);
        }
      }

      return do_put_tm(a_next, a_ios, a_fill, 
                       to_tm(a_time), format);
    }
    
    OutItrT put(OutItrT next, std::ios_base& a_ios, 
                char_type fill, const period_type& p) const 
    {
      return m_period_formatter.put_period(next, a_ios, fill,p,*this);
    }


  protected:

    static 
    string_type 
    fractional_seconds_as_string(const time_duration_type& a_time,
                                 bool null_when_zero) 
    {
      typename time_duration_type::fractional_seconds_type frac_sec = 
        a_time.fractional_seconds();

      if (null_when_zero && (frac_sec == 0)) {
        return string_type();
      }

      std::basic_ostringstream<char_type> ss;
      ss.imbue(std::locale::classic()); // don't want any formatting
      ss << std::setw(time_duration_type::num_fractional_digits())
         << std::setfill(static_cast<char_type>('0'));
#if (defined(BOOST_MSVC) && (_MSC_VER <= 1200))  // 1200 == VC++ 6.0
      // JDG [7/6/02 VC++ compatibility]
      char_type buff[34];
      ss << _i64toa(static_cast<boost::int64_t>(frac_sec), buff, 10);
#else
      ss << frac_sec;
#endif
      return ss.str();
    }

  private:
    //! Use the std facet to figure out the seperator
    /*! Simple algorithm here -- print as %X time and search for 
     *  non-numeric character.  This is assumed to be seperator. 
     *  Would be better if the standard facet was enhanced.
     */
    void set_duration_seperator_from_facet() 
    {
      //todo write this function
    }
    string_type m_time_duration_seperator;

  };
  
  template <class time_type, class CharT, class OutItrT>  
  std::locale::id time_facet<time_type, CharT, OutItrT>::id;

  template <class time_type, class CharT, class OutItrT>  
  const typename time_facet<time_type, CharT, OutItrT>::char_type 
  time_facet<time_type, CharT, OutItrT>::fractional_seconds_format[3] = {'%','f'};

  template <class time_type, class CharT, class OutItrT>  
  const typename time_facet<time_type, CharT, OutItrT>::char_type 
  time_facet<time_type, CharT, OutItrT>::fractional_seconds_or_none_format[3] = {'%','F'};

  template <class time_type, class CharT, class OutItrT>  
  const typename time_facet<time_type, CharT, OutItrT>::char_type 
  time_facet<time_type, CharT, OutItrT>::seconds_with_fractional_seconds_format[3] = 
    {'%','s'};


  template <class time_type, class CharT, class OutItrT>  
  const typename time_facet<time_type, CharT, OutItrT>::char_type 
  time_facet<time_type, CharT, OutItrT>::seconds_format[3] =  {'%','S'};

  template <class time_type, class CharT, class OutItrT>  
  const typename time_facet<time_type, CharT, OutItrT>::char_type 
  time_facet<time_type, CharT, OutItrT>::standard_format[3] =  {'%','c'};

  template <class time_type, class CharT, class OutItrT>  
  const typename time_facet<time_type, CharT, OutItrT>::char_type 
  time_facet<time_type, CharT, OutItrT>::duration_seperator[2] =  {':'};

  template <class time_type, class CharT, class OutItrT>  
  const typename time_facet<time_type,CharT, OutItrT>::char_type 
  time_facet<time_type,CharT, OutItrT>::iso_time_format_specifier[16] =  
    {'%', 'Y', '%', 'm', '%', 'd', 'T', '%', 'H', '%', 'M', '%', 'S', '%', 'F'};

  template <class time_type, class CharT, class OutItrT>  
  const typename time_facet<time_type, CharT, OutItrT>::char_type 
  time_facet<time_type, CharT, OutItrT>::iso_time_format_extended_specifier[20] =  
    {'%', 'Y', '-', '%', 'm', '-', '%', 'd', ' ', '%', 'H', ':', '%', 'M', ':', '%', 'S', '%', 'F'};


} }


#endif
