//-----------------------------------------------------------------------------
// boost mpl/aux_/preprocessor/partial_spec_params.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_PREPROCESSOR_PARTIAL_SPEC_PARAMS_HPP_INCLUDED
#define BOOST_MPL_AUX_PREPROCESSOR_PARTIAL_SPEC_PARAMS_HPP_INCLUDED

#include "boost/mpl/limits/arity.hpp"
#include "boost/mpl/aux_/preprocessor/params.hpp"
#include "boost/mpl/aux_/preprocessor/enum.hpp"
#include "boost/mpl/aux_/preprocessor/sub.hpp"
#include "boost/preprocessor/comma_if.hpp"

#define BOOST_MPL_PP_PARTIAL_SPEC_PARAMS(n, param, def) \
BOOST_MPL_PP_PARAMS(n, param) \
BOOST_PP_COMMA_IF(BOOST_MPL_PP_SUB(BOOST_MPL_METAFUNCTION_MAX_ARITY,n)) \
BOOST_MPL_PP_ENUM( \
      BOOST_MPL_PP_SUB(BOOST_MPL_METAFUNCTION_MAX_ARITY,n) \
    , def \
    ) \
/**/

#endif // BOOST_MPL_AUX_PREPROCESSOR_PARTIAL_SPEC_PARAMS_HPP_INCLUDED
