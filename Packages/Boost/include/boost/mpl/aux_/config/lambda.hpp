//-----------------------------------------------------------------------------
// boost mpl/aux_/config/lambda.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_CONFIG_LAMBDA_HPP_INCLUDED
#define BOOST_MPL_AUX_CONFIG_LAMBDA_HPP_INCLUDED

#include "boost/mpl/aux_/config/ttp.hpp"
#include "boost/config.hpp"

// agurt, 15/jan/02: full-fledged implementation requires both 
// template template parameters _and_ partial specialization

#if defined(BOOST_NO_TEMPLATE_TEMPLATE_PARAMETERS) \
 || defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) \
 && !defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT)
#   define BOOST_MPL_NO_FULL_LAMBDA_SUPPORT
#endif

//#define BOOST_MPL_NO_UNNAMED_PLACEHOLDER_SUPPORT
//#define BOOST_MPL_NO_LAMBDA_HEURISTIC

#endif // BOOST_MPL_AUX_CONFIG_LAMBDA_HPP_INCLUDED
