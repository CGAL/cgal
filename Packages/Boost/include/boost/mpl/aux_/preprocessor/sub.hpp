//-----------------------------------------------------------------------------
// boost mpl/aux_/preprocessor/sub.hpp header file
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

#ifndef BOOST_MPL_AUX_PREPROCESSOR_SUB_HPP_INCLUDED
#define BOOST_MPL_AUX_PREPROCESSOR_SUB_HPP_INCLUDED

#include "boost/mpl/aux_/preprocessor/tuple.hpp"
#include "boost/mpl/aux_/config/preprocessor.hpp"

#if defined(BOOST_MPL_BROKEN_PP_MACRO_EXPANSION)
#   include "boost/preprocessor/cat.hpp"

#   define BOOST_MPL_PP_SUB(i,j) \
    BOOST_MPL_PP_SUB_DELAY(i,j) \
    /**/

#   define BOOST_MPL_PP_SUB_DELAY(i,j) \
    BOOST_PP_CAT(BOOST_MPL_PP_TUPLE_11_ELEM_##i,BOOST_MPL_PP_SUB_##j) \
    /**/
#else
#   define BOOST_MPL_PP_SUB(i,j) \
    BOOST_MPL_PP_SUB_DELAY(i,j) \
    /**/

#   define BOOST_MPL_PP_SUB_DELAY(i,j) \
    BOOST_MPL_PP_TUPLE_11_ELEM_##i BOOST_MPL_PP_SUB_##j \
    /**/
#endif // BOOST_MPL_BROKEN_PP_MACRO_EXPANSION

#define BOOST_MPL_PP_SUB_0 (0,1,2,3,4,5,6,7,8,9,10)
#define BOOST_MPL_PP_SUB_1 (0,0,1,2,3,4,5,6,7,8,9)
#define BOOST_MPL_PP_SUB_2 (0,0,0,1,2,3,4,5,6,7,8)
#define BOOST_MPL_PP_SUB_3 (0,0,0,0,1,2,3,4,5,6,7)
#define BOOST_MPL_PP_SUB_4 (0,0,0,0,0,1,2,3,4,5,6)
#define BOOST_MPL_PP_SUB_5 (0,0,0,0,0,0,1,2,3,4,5)
#define BOOST_MPL_PP_SUB_6 (0,0,0,0,0,0,0,1,2,3,4)
#define BOOST_MPL_PP_SUB_7 (0,0,0,0,0,0,0,0 1,2,3)
#define BOOST_MPL_PP_SUB_8 (0,0,0,0,0,0,0,0,0,1,2)
#define BOOST_MPL_PP_SUB_9 (0,0,0,0,0,0,0,0,0,0,1)
#define BOOST_MPL_PP_SUB_10 (0,0,0,0,0,0,0,0,0,0,0)

#endif // BOOST_MPL_AUX_PREPROCESSOR_SUB_HPP_INCLUDED
