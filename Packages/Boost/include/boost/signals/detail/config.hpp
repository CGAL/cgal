/*
 *
 * Copyright (c) 1998-2002
 * Dr John Maddock
 *
 * Copyright (c) 2003
 * Doug Gregor
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Dr John Maddock makes no representations
 * about the suitability of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 *
 */

#ifndef BOOST_SIGNALS_CONFIG_HPP
#define BOOST_SIGNALS_CONFIG_HPP

#ifdef BOOST_HAS_DECLSPEC
#  if defined(BOOST_ALL_DYN_LINK) || defined(BOOST_SIGNALS_DYN_LINK)
#    ifdef BOOST_SIGNALS_SOURCE
#      define BOOST_SIGNALS_DECL __declspec(dllexport)
#    else
#      define BOOST_SIGNALS_DECL __declspec(dllimport)
#    endif  // BOOST_SIGNALS_SOURCE
#  endif  // DYN_LINK
#endif  // BOOST_HAS_DECLSPEC

#ifndef BOOST_SIGNALS_DECL
#  define BOOST_SIGNALS_DECL
#endif

// Setup autolinking
#if !defined(BOOST_SIGNALS_SOURCE) && !defined(BOOST_ALL_NO_LIB) && !defined(BOOST_SIGNALS_NO_LIB)
#  define BOOST_LIB_NAME boost_signals

#  if defined(BOOST_ALL_DYN_LINK) || defined(BOOST_SIGNALS_DYN_LINK)
#    define BOOST_DYN_LINK
#  endif

#  include <boost/config/auto_link.hpp>
#endif // autolinking on

#endif // BOOST_SIGNALS_CONFIG_HPP









