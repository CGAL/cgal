#ifndef _DATE_TIME_TIME_PARSING_HPP___
#define _DATE_TIME_TIME_PARSING_HPP___

/* Copyright (c) 2002,2003 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst
 * $Date$
 */

#include "boost/tokenizer.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/date_time/date_parsing.hpp"
#include "boost/cstdint.hpp"
#include <iostream>

namespace boost {
namespace date_time {

  //! Creates a time_duration object from a delimited string
  /*! Expected format for string is "[-]h[h][:mm][:ss][.fff]".
   * A negative duration will be created if the first character in
   * string is a '-', all other '-' will be treated as delimiters.
   * Accepted delimiters are "-:,.". */
  template<class time_duration>
  inline
  time_duration
  parse_delimited_time_duration(const std::string& s)
  {
    unsigned short min=0, sec =0;
    int hour =0; 
    bool is_neg = (s.at(0) == '-');
    boost::int64_t fs=0;
    int pos = 0;
    
    char_separator<char> sep("-:,.");
    tokenizer<char_separator<char> > tok(s,sep);
    for(tokenizer<char_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg){
      switch(pos) {
      case 0: {
        hour = boost::lexical_cast<int>(*beg);
        break;
      }
      case 1: {
        min = boost::lexical_cast<unsigned short>(*beg);
        break;
      }
      case 2: {
        sec = boost::lexical_cast<unsigned short>(*beg);
        break;
      };
      case 3: {
        //Works around a bug in MSVC 6 library that does not support
        //operator>> thus meaning lexical_cast will fail to compile.
#if (defined(BOOST_MSVC) && (_MSC_VER <= 1200))  // 1200 == VC++ 6.0
        fs = _atoi64(beg->c_str());
#else
        fs = boost::lexical_cast<boost::int64_t>(*beg);
#endif
        break;
      }
      }//switch
      pos++;
    }
    if(is_neg) {
      return -time_duration(hour, min, sec, fs);
    }
    else {
      return time_duration(hour, min, sec, fs);
    }
  }

  //! Utility function to split appart string
  inline
  bool 
  split(const std::string& s,
        char sep,
        std::string& first,
        std::string& second)
  {
    int sep_pos = static_cast<int>(s.find(sep));
    first = s.substr(0,sep_pos);
    second = s.substr(sep_pos+1);
    return true;
  }


  template<class time_type>
  inline
  time_type
  parse_delimited_time(const std::string& s, char sep)
  {
    typedef typename time_type::time_duration_type time_duration;
    typedef typename time_type::date_type date_type;

    //split date/time on a unique delimiter char such as ' ' or 'T'
    std::string date_string, tod_string;
    split(s, sep, date_string, tod_string);
    //call parse_date with first string
    date_type d = parse_date<date_type>(date_string);
    //call parse_time_duration with remaining string
    time_duration td = parse_delimited_time_duration<time_duration>(tod_string);
    //construct a time
    return time_type(d, td);

  }

  //! Parse time duration part of an iso time of form: [-]hhmmss (eg: 120259 is 12 hours 2 min 59 seconds)
  template<class time_duration>
  inline
  time_duration
  parse_undelimited_time_duration(const std::string& s)
  {
    int offsets[] = {2,2,2};
    int pos = 0, sign = 0;
    int hours = 0;
    short min=0, sec=0;
    // increment one position if the string was "signed"
    if(s.at(sign) == '-')
    {
      ++sign;
    }
    // stlport choked when passing s.substr() to tokenizer
    // using a new string fixed the error
    std::string remain = s.substr(sign);
    boost::offset_separator osf(offsets, offsets+3); 
    boost::tokenizer<boost::offset_separator> tok(remain, osf);
    for(boost::tokenizer<boost::offset_separator>::iterator ti=tok.begin(); ti!=tok.end();++ti){
      switch(pos) {
      case 0: 
        {
          hours = boost::lexical_cast<int>(*ti); 
          break;
        }
      case 1: 
        {
          min = boost::lexical_cast<short>(*ti); 
          break;
        }
      case 2: 
       {
         sec = boost::lexical_cast<short>(*ti); 
         break;
        }
      };
      pos++;
    }
    if(sign) {
      return -time_duration(hours, min, sec);
    }
    else {
      return time_duration(hours, min, sec);
    }
  }

  //! Parse time string of form YYYYMMDDThhmmss where T is delimeter between date and time
  template<class time_type>
  inline
  time_type
  parse_iso_time(const std::string& s, char sep)
  {
    typedef typename time_type::time_duration_type time_duration;
    typedef typename time_type::date_type date_type;

    //split date/time on a unique delimiter char such as ' ' or 'T'
    std::string date_string, tod_string;
    split(s, sep, date_string, tod_string);
    //call parse_date with first string
    date_type d = parse_undelimited_date<date_type>(date_string);
    //call parse_time_duration with remaining string
    time_duration td = parse_undelimited_time_duration<time_duration>(tod_string);
    //construct a time
    return time_type(d, td);
  }



} }//namespace date_time




#endif
