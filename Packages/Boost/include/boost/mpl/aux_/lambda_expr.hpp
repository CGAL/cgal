//-----------------------------------------------------------------------------
// boost mpl/aux_/lambda_expr.hpp header file
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

#ifndef BOOST_MPL_AUX_LAMBDA_EXPR_HPP_INCLUDED
#define BOOST_MPL_AUX_LAMBDA_EXPR_HPP_INCLUDED

#include "boost/mpl/aux_/config/lambda.hpp"

#if defined(BOOST_MPL_NO_LAMBDA_HEURISTIC)
#   define BOOST_MPL_AUX_IS_LAMBDA_EXPR(value) /**/
#else
#   define BOOST_MPL_AUX_IS_LAMBDA_EXPR(value) \
    typedef value is_le; \
    /**/
#endif

#endif // BOOST_MPL_AUX_LAMBDA_EXPR_HPP_INCLUDED
