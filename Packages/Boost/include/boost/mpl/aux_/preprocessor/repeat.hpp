//-----------------------------------------------------------------------------
// boost mpl/aux_/preprocessor/repeat.hpp header file
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

#ifndef BOOST_MPL_AUX_PREPROCESSOR_REPEAT_HPP_INCLUDED
#define BOOST_MPL_AUX_PREPROCESSOR_REPEAT_HPP_INCLUDED

#include "boost/mpl/aux_/config/preprocessor.hpp"

#if !defined(BOOST_MPL_NO_OWN_PP_PRIMITIVES)

#   include "boost/preprocessor/cat.hpp"

#   define BOOST_MPL_PP_REPEAT(n,f,param) \
    BOOST_PP_CAT(BOOST_MPL_PP_REPEAT_,n)(f,param) \
    /**/
    
#   define BOOST_MPL_PP_REPEAT_0(f,p)
#   define BOOST_MPL_PP_REPEAT_1(f,p) f(0,0,p)
#   define BOOST_MPL_PP_REPEAT_2(f,p) f(0,0,p) f(0,1,p)
#   define BOOST_MPL_PP_REPEAT_3(f,p) f(0,0,p) f(0,1,p) f(0,2,p)
#   define BOOST_MPL_PP_REPEAT_4(f,p) f(0,0,p) f(0,1,p) f(0,2,p) f(0,3,p)
#   define BOOST_MPL_PP_REPEAT_5(f,p) f(0,0,p) f(0,1,p) f(0,2,p) f(0,3,p) f(0,4,p)
#   define BOOST_MPL_PP_REPEAT_6(f,p) f(0,0,p) f(0,1,p) f(0,2,p) f(0,3,p) f(0,4,p) f(0,5,p)
#   define BOOST_MPL_PP_REPEAT_7(f,p) f(0,0,p) f(0,1,p) f(0,2,p) f(0,3,p) f(0,4,p) f(0,5,p) f(0,6,p)
#   define BOOST_MPL_PP_REPEAT_8(f,p) f(0,0,p) f(0,1,p) f(0,2,p) f(0,3,p) f(0,4,p) f(0,5,p) f(0,6,p) f(0,7,p)
#   define BOOST_MPL_PP_REPEAT_9(f,p) f(0,0,p) f(0,1,p) f(0,2,p) f(0,3,p) f(0,4,p) f(0,5,p) f(0,6,p) f(0,7,p) f(0,8,p)
#   define BOOST_MPL_PP_REPEAT_10(f,p) f(0,0,p) f(0,1,p) f(0,2,p) f(0,3,p) f(0,4,p) f(0,5,p) f(0,6,p) f(0,7,p) f(0,8,p) f(0,9,p)

#else

#   include "boost/preprocessor/repeat.hpp"

#   define BOOST_MPL_PP_REPEAT(n,f,param) \
    BOOST_PP_REPEAT_1(n,f,param) \
    /**/

#endif // BOOST_MPL_NO_OWN_PP_PRIMITIVES

#endif // BOOST_MPL_AUX_PREPROCESSOR_REPEAT_HPP_INCLUDED
