//-----------------------------------------------------------------------------
// boost mpl/aux_/arg_typedef.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2001-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_ARG_TYPEDEF_HPP_INCLUDED
#define BOOST_MPL_AUX_ARG_TYPEDEF_HPP_INCLUDED

#include "boost/mpl/aux_/config/lambda.hpp"

#if defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT)
#   define BOOST_MPL_AUX_ARG_TYPEDEF(T, name) typedef T name;
#else
#   define BOOST_MPL_AUX_ARG_TYPEDEF(T, name) /**/
#endif

#endif // BOOST_MPL_AUX_ARG_TYPEDEF_HPP_INCLUDED
