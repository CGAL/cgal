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
  *   FILE         regex_format.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Provides formatting output routines for search and replace
  *                operations.  Note this is an internal header file included
  *                by regex.hpp, do not include on its own.
  */

#ifndef BOOST_REGEX_V4_REGEX_REPLACE_HPP
#define BOOST_REGEX_V4_REGEX_REPLACE_HPP


namespace boost{

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

template <class OutputIterator, class Iterator, class traits, class Allocator, class charT>
OutputIterator regex_replace(OutputIterator out,
                         Iterator first,
                         Iterator last,
                         const reg_expression<charT, traits, Allocator>& e, 
                         const charT* fmt, 
                         match_flag_type flags = match_default)
{
   Iterator l = first;
   re_detail::merge_out_predicate<OutputIterator, Iterator, charT, Allocator, traits> oi(out, l, fmt, flags, e.get_traits());
   regex_grep(oi, first, last, e, flags);
   return (flags & format_no_copy) ? out : re_detail::re_copy_out(out, l, last);
}

template <class OutputIterator, class Iterator, class traits, class Allocator, class charT>
inline OutputIterator regex_replace(OutputIterator out,
                         Iterator first,
                         Iterator last,
                         const reg_expression<charT, traits, Allocator>& e, 
                         const std::basic_string<charT>& fmt,
                         match_flag_type flags = match_default)
{
   return regex_replace(out, first, last, e, fmt.c_str(), flags);
}

template <class traits, class Allocator, class charT>
std::basic_string<charT> regex_replace(const std::basic_string<charT>& s,
                         const reg_expression<charT, traits, Allocator>& e, 
                         const charT* fmt,
                         match_flag_type flags = match_default)
{
   std::basic_string<charT> result;
   re_detail::string_out_iterator<std::basic_string<charT> > i(result);
   regex_replace(i, s.begin(), s.end(), e, fmt, flags);
   return result;
}

template <class traits, class Allocator, class charT>
std::basic_string<charT> regex_replace(const std::basic_string<charT>& s,
                         const reg_expression<charT, traits, Allocator>& e, 
                         const std::basic_string<charT>& fmt,
                         match_flag_type flags = match_default)
{
   std::basic_string<charT> result;
   re_detail::string_out_iterator<std::basic_string<charT> > i(result);
   regex_replace(i, s.begin(), s.end(), e, fmt.c_str(), flags);
   return result;
}

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

} // namespace boost

#endif  // BOOST_REGEX_V4_REGEX_REPLACE_HPP


