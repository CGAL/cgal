#ifndef BOOST_DATE_TIME_STRING_PARSE_TREE___HPP__
#define BOOST_DATE_TIME_STRING_PARSE_TREE___HPP__

/* Copyright (c) 2004 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland
 * $Date$
 */


#include "boost/lexical_cast.hpp" //error without?
#include "boost/algorithm/string.hpp" //todo narrow this
#include <map>
#include <string>
#include <vector>
#include <algorithm>

namespace boost { namespace date_time {


template<typename charT>
struct parse_match_result
{
  parse_match_result() :
    match_depth(0),
    current_match(0)
  {}
  typedef std::basic_string<charT> string_type;
  string_type remaining() const
  {
    if (match_depth == cache.size()) {
      return string_type();
    }
    if (current_match == 0) {
      return cache;
    }
    //some of the cache was used return the rest
    return string_type(cache, match_depth); 
  }
  charT last_char() const
  {
    return cache[cache.size()-1];
  }
  bool has_remaining() const
  {
    return (cache.size() > match_depth);
  }
  string_type cache;
  unsigned short match_depth;
  unsigned short current_match;
};

  //for debug -- really only char streams...
template<typename charT>
std::basic_ostream<charT>&
operator<<(std::basic_ostream<charT>& os, parse_match_result<charT>& mr)
{
  os << "cm: " << mr.current_match 
     << " C: '" << mr.cache 
     << "' md: " << mr.match_depth
     << " R: " << mr.remaining();
  return os;
}



//! Recursive data structure to allow efficient parsing of various strings
/*! This class provides a quick lookup by building what amounts to a
 *  tree data structure.  It also features a match function which can
 *  can handle nasty input interators by caching values as it recurses
 *  the tree so that it can backtrack as needed.
 */
template<typename charT>
struct string_parse_tree
{
  typedef std::multimap<charT, string_parse_tree> ptree_coll;
  typedef typename ptree_coll::value_type value_type;
  typedef typename ptree_coll::iterator iterator;
  typedef typename ptree_coll::const_iterator const_iterator;
  typedef std::basic_string<charT> string_type;
  typedef std::vector<std::basic_string<charT> > collection_type;
  typedef parse_match_result<charT> parse_match_result_type;

  string_parse_tree(collection_type names)
  {
    // iterate thru all the elements and build the tree
    unsigned short index = 0;
    while (index != names.size() ) {
      string_type s = boost::algorithm::to_lower_copy(names[index]);
      insert(s, index+1);
      index++;
    }
    //set the last tree node = index+1  indicating a value
    index++;
  }


  string_parse_tree(unsigned short value = 0) :
    m_value(value)
  {}
  ptree_coll m_next_chars;
  unsigned short m_value;

  void insert(const string_type& s, unsigned short value)
  {
    int i = 0;
    iterator ti;
    while(i < s.size()) {
      if (i==0) {
        if (i == (s.size()-1)) {
          ti = m_next_chars.insert(value_type(s[i], 
                                              string_parse_tree<charT>(value)));
        }
        else {
          ti = m_next_chars.insert(value_type(s[i], 
                                              string_parse_tree<charT>()));
        }
      }
      else {
        if (i == (s.size()-1)) {
          ti = ti->second.m_next_chars.insert(value_type(s[i], 
                                                         string_parse_tree<charT>(value)));
        }
        
        else {
          ti = ti->second.m_next_chars.insert(value_type(s[i], 
                                                         string_parse_tree<charT>()));
        }
      
      } 
      i++;
    }
  }
  
  unsigned short
  match(std::istreambuf_iterator<charT>& sitr, 
        std::istreambuf_iterator<charT>& stream_end,
        parse_match_result_type& result,
        unsigned int& level)  const
  {
    
    unsigned int found = 0;
    level++;
    charT c;
    if (level > result.cache.size()) {
      if (sitr == stream_end) return 0; //bail stream exhausted
      c = std::tolower(*sitr);
      result.cache += c;
      sitr++;
    }
    else {
      c = result.cache[level-1];
    }
//     std::cout << "Level: " << level 
//               << " char: " << c 
//               << " result.cache: " << result.cache << std::endl;
    const_iterator litr = m_next_chars.lower_bound(c);
    const_iterator uitr = m_next_chars.upper_bound(c);
    while (litr != uitr) {
//       std::cout << "Level: " << level 
//                 << " char: " << c 
//                 << " result.cache: " << result.cache 
//                 << " value: " << litr->second.m_value 
//                 << " cm   : " << result.current_match 
//                 << " " << result.match_depth
//                 << std::endl;
      if (litr->second.m_value != 0) {
        if (result.match_depth < level) {
          result.current_match = litr->second.m_value;
          result.match_depth = level;
//           std::cout << " ....reset cm: " << result.current_match << std::endl;
        }
        found = litr->second.match(sitr, stream_end, 
                                   result, level);
        level--;
        //        found = result.current_match;
        //return litr->second.m_value; //found match
      }
      else {
        found = litr->second.match(sitr, stream_end, 
                                   result, level);
        level--;
      }
//       std::cout << "Level: " << level 
//                 << " char: " << c 
//                 << " result.cache: " << result.cache 
//                 << " found: " << found 
//                 << std::endl;
//       if (found != 0) {
//         std::cout << "Found bailing Level: " << level 
//                   << " char: " << c 
//                   << " result.cache: " << result.cache 
//                   << " found: " << found
//                   << std::endl;
//         return found; //bail we matched
//       }
      litr++;
    }
    return result.current_match;


  }


  parse_match_result_type
  match(std::istreambuf_iterator<charT>& sitr, 
        std::istreambuf_iterator<charT>& stream_end) const
  {
    // lookup to_lower of char in tree.
    unsigned int level = 0;
    //    string_type cache;
    parse_match_result_type result;
    match(sitr, stream_end, result, level);
    return result;
  }

  void printme(std::ostream& os, int& level)
  {
    level++;
    iterator itr = m_next_chars.begin();
    iterator end = m_next_chars.end();
    //    os << "starting level: " << level << std::endl;
    while (itr != end) {
      os << "level:  " << level 
         << " node:  " << itr->first 
         << " value: " << itr->second.m_value
         << std::endl;
      itr->second.printme(os, level);
      itr++;
    }
    level--;
  }

  void print(std::ostream& os)
  {
    int level = 0;
    printme(os, level);
  }
    
  void printmatch(std::ostream& os, charT c)
  {
    iterator litr = m_next_chars.lower_bound(c);
    iterator uitr = m_next_chars.upper_bound(c);
    os << "matches for: " << c << std::endl;
    while (litr != uitr) {
      os << " node:  " << litr->first 
         << " value: " << litr->second.m_value
         << std::endl;
      litr++;
    }
  }

};


} } //namespace
#endif
