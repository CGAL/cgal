//-----------------------------------------------------------------------------
// boost mpl/aux_/preprocessor/params.hpp header file
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

#ifndef BOOST_MPL_AUX_PREPROCESSOR_PARAMS_HPP_INCLUDED
#define BOOST_MPL_AUX_PREPROCESSOR_PARAMS_HPP_INCLUDED

#include "boost/mpl/aux_/config/preprocessor.hpp"

// BOOST_MPL_PP_PARAMS(0,T): <nothing>
// BOOST_MPL_PP_PARAMS(1,T): T1
// BOOST_MPL_PP_PARAMS(2,T): T1, T2
// BOOST_MPL_PP_PARAMS(n,T): T1, T2, .., Tn

#if !defined(BOOST_MPL_NO_OWN_PP_PRIMITIVES)

#   include "boost/preprocessor/cat.hpp"

#   define BOOST_MPL_PP_PARAMS(n,p) \
    BOOST_PP_CAT(BOOST_MPL_PP_PARAMS_,n)(p) \
    /**/

#   define BOOST_MPL_PP_PARAMS_0(p)
#   define BOOST_MPL_PP_PARAMS_1(p) p##1
#   define BOOST_MPL_PP_PARAMS_2(p) p##1,p##2
#   define BOOST_MPL_PP_PARAMS_3(p) p##1,p##2,p##3
#   define BOOST_MPL_PP_PARAMS_4(p) p##1,p##2,p##3,p##4
#   define BOOST_MPL_PP_PARAMS_5(p) p##1,p##2,p##3,p##4,p##5
#   define BOOST_MPL_PP_PARAMS_6(p) p##1,p##2,p##3,p##4,p##5,p##6
#   define BOOST_MPL_PP_PARAMS_7(p) p##1,p##2,p##3,p##4,p##5,p##6,p##7
#   define BOOST_MPL_PP_PARAMS_8(p) p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8
#   define BOOST_MPL_PP_PARAMS_9(p) p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9

#else

#   include "boost/preprocessor/comma_if.hpp"
#   include "boost/preprocessor/repeat.hpp"
#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/cat.hpp"

#   define BOOST_MPL_PP_AUX_PARAM_FUNC(unused, i, param) \
    BOOST_PP_COMMA_IF(i) \
    BOOST_PP_CAT(param, BOOST_PP_INC(i)) \
    /**/

#   define BOOST_MPL_PP_PARAMS(n, param) \
    BOOST_PP_REPEAT_1( \
          n \
        , BOOST_MPL_PP_AUX_PARAM_FUNC \
        , param \
        ) \
    /**/

#endif // BOOST_MPL_NO_OWN_PP_PRIMITIVES

#endif // BOOST_MPL_AUX_PREPROCESSOR_PARAMS_HPP_INCLUDED
