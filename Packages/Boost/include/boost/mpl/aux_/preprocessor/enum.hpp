//-----------------------------------------------------------------------------
// boost mpl/aux_/preprocessor/enum.hpp header file
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

#ifndef BOOST_MPL_AUX_PREPROCESSOR_ENUM_HPP_INCLUDED
#define BOOST_MPL_AUX_PREPROCESSOR_ENUM_HPP_INCLUDED

#include "boost/mpl/aux_/config/preprocessor.hpp"

// BOOST_MPL_PP_ENUM(0,int): <nothing>
// BOOST_MPL_PP_ENUM(1,int): int
// BOOST_MPL_PP_ENUM(2,int): int, int
// BOOST_MPL_PP_ENUM(n,int): int, int, .., int

#if !defined(BOOST_MPL_NO_OWN_PP_PRIMITIVES)

#   include "boost/preprocessor/cat.hpp"

#   define BOOST_MPL_PP_ENUM(n, param) \
    BOOST_PP_CAT(BOOST_MPL_PP_ENUM_,n)(param) \
    /**/
    
#   define BOOST_MPL_PP_ENUM_0(p)
#   define BOOST_MPL_PP_ENUM_1(p) p
#   define BOOST_MPL_PP_ENUM_2(p) p,p
#   define BOOST_MPL_PP_ENUM_3(p) p,p,p
#   define BOOST_MPL_PP_ENUM_4(p) p,p,p,p
#   define BOOST_MPL_PP_ENUM_5(p) p,p,p,p,p
#   define BOOST_MPL_PP_ENUM_6(p) p,p,p,p,p,p
#   define BOOST_MPL_PP_ENUM_7(p) p,p,p,p,p,p,p
#   define BOOST_MPL_PP_ENUM_8(p) p,p,p,p,p,p,p,p
#   define BOOST_MPL_PP_ENUM_9(p) p,p,p,p,p,p,p,p,p

#else

#   include "boost/preprocessor/comma_if.hpp"
#   include "boost/preprocessor/repeat.hpp"

#   define BOOST_MPL_PP_AUX_ENUM_FUNC(unused, i, param) \
    BOOST_PP_COMMA_IF(i) param \
    /**/

#   define BOOST_MPL_PP_ENUM(n, param) \
    BOOST_PP_REPEAT_1( \
          n \
        , BOOST_MPL_PP_AUX_ENUM_FUNC \
        , param \
        ) \
    /**/

#endif // BOOST_MPL_NO_OWN_PP_PRIMITIVES

#endif // BOOST_MPL_AUX_PREPROCESSOR_ENUM_HPP_INCLUDED
