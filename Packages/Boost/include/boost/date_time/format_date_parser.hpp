
#ifndef DATE_TIME_FORMAT_DATE_PARSER_HPP__
#define DATE_TIME_FORMAT_DATE_PARSER_HPP__

/* Copyright (c) 2004 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland
 * $Date$
 */


#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp" //todo narrow this
#include "boost/date_time/string_parse_tree.hpp"
#include "boost/date_time/strings_from_facet.hpp"
#include <string>
#include <vector>

namespace boost { namespace date_time {

//!Helper function for parsing fixed length strings into integers
template<typename int_type, typename charT>
inline
int_type
fixed_string_to_int(std::istreambuf_iterator<charT>& itr,
                    unsigned int length)
{
  typedef std::basic_string<charT>  string_type;
  unsigned int j = 0;
  string_type s;
  while (j < length) {
    s += (*itr);
    itr++;
    j++;
  }
  int_type i = boost::lexical_cast<int_type>(s);
  return i;
}

//!Helper function for parsing fixed length strings into integers
template<typename int_type, typename charT>
inline
int_type
var_string_to_int(std::istreambuf_iterator<charT>& itr,
                  unsigned int max_length)
{
  typedef std::basic_string<charT>  string_type;
  unsigned int j = 0;
  string_type s;
  while ((j < max_length) && std::isdigit(*itr)) {
    s += (*itr);
    itr++;
    j++;
  }
  int_type i = boost::lexical_cast<int_type>(s);
  return i;
}


//! Class with generic date parsing using a format string
/*! The following is the set of recognized format specifiers
 -  %a - Short weekday name
 -  %A - Long weekday name
 -  %b - Abbreviated month name
 -  %B - Full month name
 -  %d - Day of the month as decimal
 -  %j - Day of year as decimal from 1 to 366
 -  %m - Month name as a decimal 01 to 12
 -  %U - Week number 00 to 53 with first Sunday as the first day of week 1?
 -  %w - Weekday as decimal number 0 to 6 where Sunday == 0
 -  %W - Week number 00 to 53 where Monday is first day of week 1
 -  %x - facet default date representation
 -  %y - Year without the century - eg: 04 for 2004
 -  %Y - Year with century 

 The weekday specifiers (%a and %A) do not add to the date construction,
 but they provide a way to skip over the weekday names for formats that
 provide them.

 todo -- Another interesting feature that this approach could provide is
         an option to fill in any missing fields with the current values
         from the clock.  So if you have %m-%d the parser would detect
         the missing year value and fill it in using the clock. 

 todo -- What to do with the %x.  %x in the classic facet is just bad...

 */
template<class date_type, typename charT>
class format_date_parser
{
 public:
  typedef std::basic_string<charT>        string_type;
  typedef std::basic_stringstream<charT>  stringstream_type;
  typedef std::istreambuf_iterator<charT> stream_itr_type;
  typedef typename string_type::const_iterator const_itr;
  typedef typename date_type::year_type  year_type;
  typedef typename date_type::month_type month_type;
  typedef typename date_type::duration_type duration_type;
  typedef string_parse_tree<charT> parse_tree_type;
  typedef typename parse_tree_type::parse_match_result_type match_results;
  typedef std::vector<std::basic_string<charT> > input_collection_type;
  
  format_date_parser(const string_type& format,
                     const input_collection_type& month_short_names,
                     const input_collection_type& month_long_names,
                     const input_collection_type& weekday_short_names,
                     const input_collection_type& weekday_long_names) :
    m_format(format),
    m_month_short_names(month_short_names),
    m_month_long_names(month_long_names),
    m_weekday_short_names(weekday_short_names),
    m_weekday_long_names(weekday_long_names)
  {}
  
  format_date_parser(const string_type& format,
                     const std::locale& locale) :
    m_format(format),
    m_month_short_names(gather_month_strings<charT>(locale)),
    m_month_long_names(gather_month_strings<charT>(locale, false)),
    m_weekday_short_names(gather_weekday_strings<char>(locale)),
    m_weekday_long_names(gather_weekday_strings<char>(locale, false))
  {}
  
  string_type format() const
  {
    return m_format;
  }

  void format(string_type format)
  {
    m_format = format;
  }

  date_type
  parse_date(const string_type& value, 
             const string_type& format) const
  {
    stringstream_type ss;
    ss << value; 
    stream_itr_type sitr(ss);
    stream_itr_type stream_end;
    return parse_date(sitr, stream_end, format);
  }

  date_type
  parse_date(std::istreambuf_iterator<charT>& sitr, 
             std::istreambuf_iterator<charT>& stream_end) const
  {
    parse_date(sitr, stream_end, m_format);
  }

  date_type
  parse_date(std::istreambuf_iterator<charT>& sitr, 
             std::istreambuf_iterator<charT>& stream_end,
             string_type format) const
  {
    bool use_current_char = false;
    charT current_char = *sitr;

    unsigned short year(0), month(0), day(0), day_of_year(0);
    
    const_itr itr(format.begin());
    while (itr != format.end() && (sitr != stream_end)) {
      if (*itr == '%') {
        itr++;
        if (*itr != '%') { //ignore '%%'
          unsigned short i = 0;
          switch(*itr) {
          case 'a': 
            {
              //this value is just throw away.  It could be used for
              //error checking potentially, but it isn't helpful in 
              //actually constructing the date - we just need to get it
              //out of the stream
              match_results mr = m_weekday_short_names.match(sitr, stream_end);
              unsigned int wkday = mr.current_match;
              if (mr.has_remaining()) {
                current_char = mr.last_char();
                use_current_char = true;
              }
              break;
            }
          case 'A': 
            {
              //this value is just throw away.  It could be used for
              //error checking potentially, but it isn't helpful in 
              //actually constructing the date - we just need to get it
              //out of the stream
              match_results mr = m_weekday_long_names.match(sitr, stream_end);
              unsigned int wkday = mr.current_match;
              if (mr.has_remaining()) {
                current_char = mr.last_char();
                use_current_char = true;
              }
              break;
            }
          case 'b': 
            {
              match_results mr = m_month_short_names.match(sitr, stream_end);
              month = mr.current_match;
              if (mr.has_remaining()) {
                current_char = mr.last_char();
                use_current_char = true;
              }
              break;
            }
          case 'B': 
            {
              match_results mr = m_month_long_names.match(sitr, stream_end);
              month = mr.current_match;
              if (mr.has_remaining()) {
                current_char = mr.last_char();
                use_current_char = true;
              }
              break;
            }
          case 'd': 
            {
              day = var_string_to_int<unsigned short, charT>(sitr, 2);
              break;
            }
          case 'j': 
            {
              day_of_year = fixed_string_to_int<unsigned short, charT>(sitr, 3);
              break;
            }
          case 'm': 
            {
              month = var_string_to_int<unsigned short, charT>(sitr, 2);
              break;
            }
          case 'Y': 
            {
              year = fixed_string_to_int<unsigned short, charT>(sitr, 4);
              break;
            }
          case 'y': 
            {
              year = fixed_string_to_int<unsigned short, charT>(sitr, 2);
              year += 2000; //make 2 digit years in this century
              break;
            }
          default:
            {} //ignore those we don't understand
            
          }//switch
          
        }
        
        itr++; //advance past format specifier
      }
      else {  //skip past chars in format and in buffer
        itr++;
        if (use_current_char) {
          use_current_char = false;
          current_char = *sitr;
        }
        else {
          sitr++;
        }
      }
    }
    
    if (day_of_year != 0) {
      date_type d(year-1,12,31); //end of prior year
      return d + duration_type(day_of_year);
    }
    return date_type(year, month, day);
  }
  
  
 private:
  string_type m_format;
  parse_tree_type m_month_short_names;
  parse_tree_type m_month_long_names;
  parse_tree_type m_weekday_short_names;
  parse_tree_type m_weekday_long_names;
};

} } //namespace

#endif
