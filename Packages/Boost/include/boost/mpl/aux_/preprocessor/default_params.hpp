//-----------------------------------------------------------------------------
// boost mpl/aux_/preprocessor/default_params.hpp header file
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

#ifndef BOOST_MPL_AUX_PREPROCESSOR_DEFAULT_PARAMS_HPP_INCLUDED
#define BOOST_MPL_AUX_PREPROCESSOR_DEFAULT_PARAMS_HPP_INCLUDED

#include "boost/mpl/aux_/config/preprocessor.hpp"

// BOOST_MPL_PP_DEFAULT_PARAMS(0,T,int): <nothing>
// BOOST_MPL_PP_DEFAULT_PARAMS(1,T,int): T1 = int
// BOOST_MPL_PP_DEFAULT_PARAMS(2,T,int): T1 = int, T2 = int
// BOOST_MPL_PP_DEFAULT_PARAMS(n,T,int): T1 = int, T2 = int, .., Tn = int

#if !defined(BOOST_MPL_NO_OWN_PP_PRIMITIVES)

#   include "boost/preprocessor/cat.hpp"

#   define BOOST_MPL_PP_DEFAULT_PARAMS(n,p,v) \
    BOOST_PP_CAT(BOOST_MPL_PP_DEFAULT_PARAMS_,n)(p,v) \
    /**/
    
#   define BOOST_MPL_PP_DEFAULT_PARAMS_0(p,v)
#   define BOOST_MPL_PP_DEFAULT_PARAMS_1(p,v) p##1=v
#   define BOOST_MPL_PP_DEFAULT_PARAMS_2(p,v) p##1=v,p##2=v
#   define BOOST_MPL_PP_DEFAULT_PARAMS_3(p,v) p##1=v,p##2=v,p##3=v
#   define BOOST_MPL_PP_DEFAULT_PARAMS_4(p,v) p##1=v,p##2=v,p##3=v,p##4=v
#   define BOOST_MPL_PP_DEFAULT_PARAMS_5(p,v) p##1=v,p##2=v,p##3=v,p##4=v,p##5=v
#   define BOOST_MPL_PP_DEFAULT_PARAMS_6(p,v) p##1=v,p##2=v,p##3=v,p##4=v,p##5=v,p##6=v
#   define BOOST_MPL_PP_DEFAULT_PARAMS_7(p,v) p##1=v,p##2=v,p##3=v,p##4=v,p##5=v,p##6=v,p##7=v
#   define BOOST_MPL_PP_DEFAULT_PARAMS_8(p,v) p##1=v,p##2=v,p##3=v,p##4=v,p##5=v,p##6=v,p##7=v,p##8=v
#   define BOOST_MPL_PP_DEFAULT_PARAMS_9(p,v) p##1=v,p##2=v,p##3=v,p##4=v,p##5=v,p##6=v,p##7=v,p##8=v,p##9=v

#else

#   include "boost/preprocessor/tuple/elem.hpp"
#   include "boost/preprocessor/comma_if.hpp"
#   include "boost/preprocessor/repeat.hpp"
#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/cat.hpp"

#   define BOOST_MPL_PP_AUX_DEFAULT_PARAM_FUNC(unused, i, pv) \
    BOOST_PP_COMMA_IF(i) \
    BOOST_PP_CAT( BOOST_PP_TUPLE_ELEM(2,0,pv), BOOST_PP_INC(i) ) \
        = BOOST_PP_TUPLE_ELEM(2,1,pv) \
    /**/

#   define BOOST_MPL_PP_DEFAULT_PARAMS(n, param, value) \
    BOOST_PP_REPEAT_1( \
          n \
        , BOOST_MPL_PP_AUX_DEFAULT_PARAM_FUNC \
        , (param,value) \
        ) \
    /**/

#endif // BOOST_MPL_USE_OWN_PP_PRIMITIVES

#endif // BOOST_MPL_AUX_PREPROCESSOR_DEFAULT_PARAMS_HPP_INCLUDED
