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
  *   FILE         regex.cpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Declares boost::reg_expression<> and associated
  *                functions and classes. This header is the main
  *                entry point for the template regex code.
  */

#ifndef BOOST_RE_REGEX_HPP_INCLUDED
#define BOOST_RE_REGEX_HPP_INCLUDED

#ifndef BOOST_RE_CREGEX_HPP
#include <boost/cregex.hpp>
#endif

#ifdef __cplusplus

// what follows is all C++ don't include in C builds!!

#ifdef BOOST_REGEX_DEBUG
# include <iosfwd>
#endif

#include <new>
#include <cstring>
#ifndef BOOST_REGEX_CONFIG_HPP
#include <boost/regex/config.hpp>
#endif
#ifndef BOOST_REGEX_FWD_HPP
#include <boost/regex_fwd.hpp>
#endif
#ifndef BOOST_REGEX_STACK_HPP
#include <boost/regex/v4/regex_stack.hpp>
#endif
#ifndef BOOST_REGEX_RAW_BUFFER_HPP
#include <boost/regex/v4/regex_raw_buffer.hpp>
#endif
#ifndef BOOST_REGEX_KMP_HPP
#include <boost/regex/v4/regex_kmp.hpp>
#endif
#ifndef BOOST_RE_PAT_EXCEPT_HPP
#include <boost/regex/pattern_except.hpp>
#endif
#ifndef BOOST_REGEX_TRAITS_HPP
#include <boost/regex/regex_traits.hpp>
#endif
#include <boost/scoped_array.hpp>

#ifndef BOOST_REGEX_V4_CHAR_REGEX_TRAITS_HPP
#include <boost/regex/v4/char_regex_traits.hpp>
#endif
#ifndef BOOST_REGEX_V4_STATES_HPP
#include <boost/regex/v4/states.hpp>
#endif
#ifndef BOOST_REGEX_V4_REGBASE_HPP
#include <boost/regex/v4/regbase.hpp>
#endif
#ifndef BOOST_REGEX_V4_ITERATOR_TRAITS_HPP
#include <boost/regex/v4/iterator_traits.hpp>
#endif
#ifndef BOOST_REGEX_V4_ITERATOR_TRAITS_HPP
#include <boost/regex/v4/iterator_traits.hpp>
#endif
#ifndef BOOST_REGEX_V4_BASIC_REGEX_HPP
#include <boost/regex/v4/basic_regex.hpp>
#endif
#ifndef BOOST_REGEX_V4_SUB_MATCH_HPP
#include <boost/regex/v4/sub_match.hpp>
#endif
#ifndef BOOST_REGEX_FORMAT_HPP
#include <boost/regex/v4/regex_format.hpp>
#endif
#ifndef BOOST_REGEX_V4_MATCH_RESULTS_HPP
#include <boost/regex/v4/match_results.hpp>
#endif
#ifndef BOOST_REGEX_COMPILE_HPP
#include <boost/regex/v4/regex_compile.hpp>
#endif

//
// template instances:
//
#define BOOST_REGEX_CHAR_T char
#ifdef BOOST_REGEX_NARROW_INSTANTIATE
#  define BOOST_REGEX_INSTANTIATE
#endif
#include <boost/regex/v4/instances.hpp>
#undef BOOST_REGEX_CHAR_T
#ifdef BOOST_REGEX_INSTANTIATE
#  undef BOOST_REGEX_INSTANTIATE
#endif

#ifndef BOOST_NO_WREGEX
#define BOOST_REGEX_CHAR_T boost::regex_wchar_type
#ifdef BOOST_REGEX_WIDE_INSTANTIATE
#  define BOOST_REGEX_INSTANTIATE
#endif
#include <boost/regex/v4/instances.hpp>
#undef BOOST_REGEX_CHAR_T
#ifdef BOOST_REGEX_INSTANTIATE
#  undef BOOST_REGEX_INSTANTIATE
#endif
#endif


namespace boost{
#ifdef BOOST_REGEX_NO_FWD
typedef reg_expression<char, regex_traits<char>, BOOST_DEFAULT_ALLOCATOR(char)> regex;
#ifndef BOOST_NO_WREGEX
typedef reg_expression<wchar_t, regex_traits<wchar_t>, BOOST_DEFAULT_ALLOCATOR(wchar_t)> wregex;
#endif
#endif

typedef match_results<const char*> cmatch;
typedef match_results<std::string::const_iterator> smatch;
#ifndef BOOST_NO_WREGEX
typedef match_results<const wchar_t*> wcmatch;
typedef match_results<std::wstring::const_iterator> wsmatch;
#endif

} // namespace boost
#ifndef BOOST_REGEX_MATCH_HPP
#include <boost/regex/v4/regex_match.hpp>
#endif
#ifndef BOOST_REGEX_V4_REGEX_SEARCH_HPP
#include <boost/regex/v4/regex_search.hpp>
#endif
#ifndef BOOST_REGEX_V4_REGEX_GREP_HPP
#include <boost/regex/v4/regex_grep.hpp>
#endif
#ifndef BOOST_REGEX_V4_REGEX_REPLACE_HPP
#include <boost/regex/v4/regex_replace.hpp>
#endif
#ifndef BOOST_REGEX_V4_REGEX_MERGE_HPP
#include <boost/regex/v4/regex_merge.hpp>
#endif
#ifndef BOOST_REGEX_SPLIT_HPP
#include <boost/regex/v4/regex_split.hpp>
#endif
#ifndef BOOST_REGEX_ITERATOR_HPP
#include <boost/regex/v4/regex_iterator.hpp>
#endif
#ifndef BOOST_REGEX_TOKEN_ITERATOR_HPP
#include <boost/regex/v4/regex_token_iterator.hpp>
#endif

#endif  // __cplusplus

#endif  // include































