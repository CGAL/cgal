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
  *   FILE         src.cpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Includes all the regex source files, include this
  *                file only if you need to build the regex library
  *                as a single file.  You must include this file
  *                before any other regex header.
  *
  *                CAUTION: THIS FILE IS DEPRICATED AND WILL CAUSE 
  *                UNNECESSARY CODE BLOAT.
  */

#if (!defined(BOOST_REGEX_NO_LIB) || !defined(BOOST_REGEX_NO_EXTERNAL_TEMPLATES)) && defined(BOOST_REGEX_CONFIG_HPP)
#error too late you have already included a regex header - make sure that you include this header before any other boost header
#endif

#define BOOST_REGEX_NO_LIB
#define BOOST_REGEX_STATIC_LINK
#define BOOST_REGEX_NO_EXTERNAL_TEMPLATES

#include <boost/regex.hpp>

//
// include library source files:
//
#ifdef BOOST_REGEX_USE_WIN32_LOCALE
#include "libs/regex/src/w32_regex_traits.cpp"
#elif defined(BOOST_REGEX_USE_C_LOCALE)
#include "libs/regex/src/c_regex_traits.cpp"
#else
#include "libs/regex/src/cpp_regex_traits.cpp"
#endif
#include "libs/regex/src/c_regex_traits_common.cpp"
#include "libs/regex/src/cregex.cpp"
#include "libs/regex/src/fileiter.cpp"
#include "libs/regex/src/posix_api.cpp"
#include "libs/regex/src/wide_posix_api.cpp"
#include "libs/regex/src/regex.cpp"
#include "libs/regex/src/regex_debug.cpp"
#include "libs/regex/src/regex_synch.cpp"



