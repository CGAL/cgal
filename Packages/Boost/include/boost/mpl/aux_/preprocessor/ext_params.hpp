//-----------------------------------------------------------------------------
// boost mpl/aux_/preprocessor/ext_params.hpp header file
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

#ifndef BOOST_MPL_AUX_PREPROCESSOR_EXT_PARAMS_HPP_INCLUDED
#define BOOST_MPL_AUX_PREPROCESSOR_EXT_PARAMS_HPP_INCLUDED

#include "boost/mpl/aux_/config/preprocessor.hpp"

// BOOST_MPL_PP_EXT_PARAMS(1,1,T): <nothing>
// BOOST_MPL_PP_EXT_PARAMS(1,2,T): T1
// BOOST_MPL_PP_EXT_PARAMS(1,3,T): T1, T2
// BOOST_MPL_PP_EXT_PARAMS(1,n,T): T1, T2, .., Tn-1

#if !defined(BOOST_MPL_NO_OWN_PP_PRIMITIVES)

#   include "boost/mpl/aux_/preprocessor/filter_params.hpp"
#   include "boost/mpl/aux_/preprocessor/sub.hpp"

#   define BOOST_MPL_PP_EXT_PARAMS(i,j,p) \
    BOOST_MPL_PP_EXT_PARAMS_DELAY_1(i,BOOST_MPL_PP_SUB(j,i),p) \
    /**/

#   define BOOST_MPL_PP_EXT_PARAMS_DELAY_1(i,n,p) \
    BOOST_MPL_PP_EXT_PARAMS_DELAY_2(i,n,p) \
    /**/

#   define BOOST_MPL_PP_EXT_PARAMS_DELAY_2(i,n,p) \
    BOOST_MPL_PP_EXT_PARAMS_##i(n,p) \
    /**/

#   define BOOST_MPL_PP_EXT_PARAMS_1(i,p) BOOST_MPL_PP_FILTER_PARAMS_##i(p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9)
#   define BOOST_MPL_PP_EXT_PARAMS_2(i,p) BOOST_MPL_PP_FILTER_PARAMS_##i(p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9,p1)
#   define BOOST_MPL_PP_EXT_PARAMS_3(i,p) BOOST_MPL_PP_FILTER_PARAMS_##i(p##3,p##4,p##5,p##6,p##7,p##8,p##9,p1,p2)
#   define BOOST_MPL_PP_EXT_PARAMS_4(i,p) BOOST_MPL_PP_FILTER_PARAMS_##i(p##4,p##5,p##6,p##7,p##8,p##9,p1,p2,p3)
#   define BOOST_MPL_PP_EXT_PARAMS_5(i,p) BOOST_MPL_PP_FILTER_PARAMS_##i(p##5,p##6,p##7,p##8,p##9,p1,p2,p3,p4)
#   define BOOST_MPL_PP_EXT_PARAMS_6(i,p) BOOST_MPL_PP_FILTER_PARAMS_##i(p##6,p##7,p##8,p##9,p1,p2,p3,p4,p5)
#   define BOOST_MPL_PP_EXT_PARAMS_7(i,p) BOOST_MPL_PP_FILTER_PARAMS_##i(p##7,p##8,p##9,p1,p2,p3,p4,p5,p6)
#   define BOOST_MPL_PP_EXT_PARAMS_8(i,p) BOOST_MPL_PP_FILTER_PARAMS_##i(p##8,p##9,p1,p2,p3,p4,p5,p6,p7)
#   define BOOST_MPL_PP_EXT_PARAMS_9(i,p) BOOST_MPL_PP_FILTER_PARAMS_##i(p##9,p1,p2,p3,p4,p5,p6,p7,p8)

#else

#   include "boost/preprocessor/arithmetic/add.hpp"
#   include "boost/preprocessor/arithmetic/sub.hpp"
#   include "boost/preprocessor/comma_if.hpp"
#   include "boost/preprocessor/repeat.hpp"
#   include "boost/preprocessor/tuple/elem.hpp"
#   include "boost/preprocessor/cat.hpp"

#   define BOOST_MPL_PP_AUX_EXT_PARAM_FUNC(unused, i, op) \
    BOOST_PP_COMMA_IF(i) \
    BOOST_PP_CAT( \
          BOOST_PP_TUPLE_ELEM(2,1,op) \
        , BOOST_PP_ADD_D(1, i, BOOST_PP_TUPLE_ELEM(2,0,op)) \
        ) \
    /**/

#   define BOOST_MPL_PP_EXT_PARAMS(i, j, param) \
    BOOST_PP_REPEAT_1( \
          BOOST_PP_SUB_D(1,j,i) \
        , BOOST_MPL_PP_AUX_EXT_PARAM_FUNC \
        , (i,param) \
        ) \
    /**/

#endif // BOOST_MPL_NO_OWN_PP_PRIMITIVES

#endif // BOOST_MPL_AUX_PREPROCESSOR_EXT_PARAMS_HPP_INCLUDED
