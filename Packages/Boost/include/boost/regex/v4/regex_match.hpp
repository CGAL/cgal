/*
 *
 * Copyright (c) 1998-2002
 * Dr John Maddock
 *
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         regex_match.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Regular expression matching algorithms.
  *                Note this is an internal header file included
  *                by regex.hpp, do not include on its own.
  */


#ifndef BOOST_REGEX_MATCH_HPP
#define BOOST_REGEX_MATCH_HPP

#ifndef BOOST_REGEX_MAX_STATE_COUNT
#  define BOOST_REGEX_MAX_STATE_COUNT 100000000
#endif

#include <boost/limits.hpp>
#include <boost/regex/v4/perl_matcher.hpp>


namespace boost{

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

//
// proc regex_match
// returns true if the specified regular expression matches
// the whole of the input.  Fills in what matched in m.
//
template <class BidiIterator, class Allocator, class charT, class traits, class Allocator2>
bool regex_match(BidiIterator first, BidiIterator last, 
                 match_results<BidiIterator, Allocator>& m, 
                 const reg_expression<charT, traits, Allocator2>& e, 
                 match_flag_type flags = match_default)
{
   re_detail::perl_matcher<BidiIterator, Allocator, traits, Allocator2> matcher(first, last, m, e, flags);
   return matcher.match();
}
template <class iterator, class charT, class traits, class Allocator2>
bool regex_match(iterator first, iterator last, 
                 const reg_expression<charT, traits, Allocator2>& e, 
                 match_flag_type flags = match_default)
{
   match_results<iterator> m;
   return regex_match(first, last, m, e, flags);
}
//
// query_match convenience interfaces:
#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
//
// this isn't really a partial specialisation, but template function
// overloading - if the compiler doesn't support partial specialisation
// then it really won't support this either:
template <class charT, class Allocator, class traits, class Allocator2>
inline bool regex_match(const charT* str, 
                        match_results<const charT*, Allocator>& m, 
                        const reg_expression<charT, traits, Allocator2>& e, 
                        match_flag_type flags = match_default)
{
   return regex_match(str, str + traits::length(str), m, e, flags);
}

template <class ST, class SA, class Allocator, class charT, class traits, class Allocator2>
inline bool regex_match(const std::basic_string<charT, ST, SA>& s, 
                 match_results<typename std::basic_string<charT, ST, SA>::const_iterator, Allocator>& m, 
                 const reg_expression<charT, traits, Allocator2>& e, 
                 match_flag_type flags = match_default)
{
   return regex_match(s.begin(), s.end(), m, e, flags);
}
template <class charT, class traits, class Allocator2>
inline bool regex_match(const charT* str, 
                        const reg_expression<charT, traits, Allocator2>& e, 
                        match_flag_type flags = match_default)
{
   match_results<const charT*> m;
   return regex_match(str, str + traits::length(str), m, e, flags);
}

template <class ST, class SA, class charT, class traits, class Allocator2>
inline bool regex_match(const std::basic_string<charT, ST, SA>& s, 
                 const reg_expression<charT, traits, Allocator2>& e, 
                 match_flag_type flags = match_default)
{
   typedef typename std::basic_string<charT, ST, SA>::const_iterator iterator;
   match_results<iterator> m;
   return regex_match(s.begin(), s.end(), m, e, flags);
}
#else  // partial ordering
inline bool regex_match(const char* str, 
                        cmatch& m, 
                        const regex& e, 
                        match_flag_type flags = match_default)
{
   return regex_match(str, str + regex::traits_type::length(str), m, e, flags);
}
inline bool regex_match(const char* str, 
                        const regex& e, 
                        match_flag_type flags = match_default)
{
   match_results<const char*> m;
   return regex_match(str, str + regex::traits_type::length(str), m, e, flags);
}
#ifndef BOOST_NO_WREGEX
inline bool regex_match(const wchar_t* str, 
                        wcmatch& m, 
                        const wregex& e, 
                        match_flag_type flags = match_default)
{
   return regex_match(str, str + wregex::traits_type::length(str), m, e, flags);
}
inline bool regex_match(const wchar_t* str, 
                        const wregex& e, 
                        match_flag_type flags = match_default)
{
   match_results<const wchar_t*> m;
   return regex_match(str, str + wregex::traits_type::length(str), m, e, flags);
}
#endif
inline bool regex_match(const std::string& s, 
                        smatch& m,
                        const regex& e, 
                        match_flag_type flags = match_default)
{
   return regex_match(s.begin(), s.end(), m, e, flags);
}
inline bool regex_match(const std::string& s, 
                        const regex& e, 
                        match_flag_type flags = match_default)
{
   match_results<std::string::const_iterator> m;
   return regex_match(s.begin(), s.end(), m, e, flags);
}
#if !defined(BOOST_NO_WREGEX)
inline bool regex_match(const std::basic_string<wchar_t>& s, 
                        wsmatch& m,
                        const wregex& e, 
                        match_flag_type flags = match_default)
{
   return regex_match(s.begin(), s.end(), m, e, flags);
}
inline bool regex_match(const std::basic_string<wchar_t>& s, 
                        const wregex& e, 
                        match_flag_type flags = match_default)
{
   match_results<std::basic_string<wchar_t>::const_iterator> m;
   return regex_match(s.begin(), s.end(), m, e, flags);
}
#endif

#endif


#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

} // namespace boost

#endif   // BOOST_REGEX_MATCH_HPP


















