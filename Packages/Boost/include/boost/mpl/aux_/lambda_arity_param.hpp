//-----------------------------------------------------------------------------
// boost mpl/aux_/lambda_arity_param.hpp header file
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

#ifndef BOOST_MPL_AUX_LAMBDA_ARITY_PARAM_HPP_INCLUDED
#define BOOST_MPL_AUX_LAMBDA_ARITY_PARAM_HPP_INCLUDED

#include "boost/mpl/aux_/config/ttp.hpp"

#if !defined(BOOST_EXTENDED_TEMPLATE_PARAMETERS_MATCHING)
#   define BOOST_MPL_AUX_LAMBDA_ARITY_PARAM(param)    
#else
#   define BOOST_MPL_AUX_LAMBDA_ARITY_PARAM(param) , param
#endif

#endif // BOOST_MPL_AUX_LAMBDA_ARITY_PARAM_HPP_INCLUDED
