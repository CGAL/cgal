//-----------------------------------------------------------------------------
// boost variant/detail/config.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2003
// Eric Friedman
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_VARIANT_DETAIL_CONFIG_HPP
#define BOOST_VARIANT_DETAIL_CONFIG_HPP

#include "boost/config.hpp"
#include "boost/detail/workaround.hpp"

///////////////////////////////////////////////////////////////////////////////
// macro BOOST_VARIANT_AUX_BROKEN_CONSTRUCTOR_TEMPLATE_ORDERING
//
#if BOOST_WORKAROUND(__MWERKS__, <= 0x3201) \
 || BOOST_WORKAROUND(BOOST_INTEL, <= 700) \
 || BOOST_WORKAROUND(BOOST_MSVC, <= 1300) \
 && !defined(BOOST_VARIANT_AUX_BROKEN_CONSTRUCTOR_TEMPLATE_ORDERING)
#   define BOOST_VARIANT_AUX_BROKEN_CONSTRUCTOR_TEMPLATE_ORDERING
#endif

///////////////////////////////////////////////////////////////////////////////
// macro BOOST_VARIANT_AUX_HAS_CONSTRUCTOR_TEMPLATE_ORDERING_SFINAE_WKND
//
#if !defined(BOOST_NO_SFINAE) \
 && !BOOST_WORKAROUND(BOOST_INTEL, <= 700) \
 && !defined(BOOST_VARIANT_AUX_HAS_CONSTRUCTOR_TEMPLATE_ORDERING_SFINAE_WKND)
#   define BOOST_VARIANT_AUX_HAS_CONSTRUCTOR_TEMPLATE_ORDERING_SFINAE_WKND
#endif

#endif // BOOST_VARIANT_DETAIL_CONFIG_HPP
