
#ifndef BOOST_MPL_AUX_PREPROCESSOR_TOKEN_EQUAL_HPP_INCLUDED
#define BOOST_MPL_AUX_PREPROCESSOR_TOKEN_EQUAL_HPP_INCLUDED

// + file: boost/mpl/aux_/preprocessor/token_equal.hpp
// + last modified: 03/may/03

// Copyright (c) 2003
// Paul Mensonides, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.
//
// See http://www.boost.org/libs/mpl for documentation.

#include "boost/mpl/aux_/preprocessor/is_seq.hpp"

#include "boost/preprocessor/if.hpp"
#include "boost/preprocessor/logical/bitand.hpp"
#include "boost/preprocessor/logical/compl.hpp"
#include "boost/preprocessor/tuple/eat.hpp"
#include "boost/preprocessor/cat.hpp"

// compares tokens 'a' and 'b' for equality:
//
//   #define BOOST_MPL_PP_TOKEN_EQUAL_apple(x) x
//   #define BOOST_MPL_PP_TOKEN_EQUAL_orange(x) x
//
//   BOOST_PP_ASSERT( BOOST_PP_NOT( BOOST_MPL_PP_TOKEN_EQUAL(apple, abc) ) )
//   BOOST_PP_ASSERT( BOOST_PP_NOT( BOOST_MPL_PP_TOKEN_EQUAL(abc, apple) ) )
//   BOOST_PP_ASSERT( BOOST_PP_NOT( BOOST_MPL_PP_TOKEN_EQUAL(apple, orange) ) )
//   BOOST_PP_ASSERT( BOOST_MPL_PP_TOKEN_EQUAL(apple, apple) )
//   BOOST_PP_ASSERT( BOOST_MPL_PP_TOKEN_EQUAL(orange, orange) )

#define BOOST_MPL_PP_TOKEN_EQUAL(a, b) \
    BOOST_PP_IIF( \
        BOOST_PP_BITAND( \
              BOOST_MPL_PP_IS_SEQ( BOOST_PP_CAT(BOOST_MPL_PP_TOKEN_EQUAL_, a)((unused)) ) \
            , BOOST_MPL_PP_IS_SEQ( BOOST_PP_CAT(BOOST_MPL_PP_TOKEN_EQUAL_, b)((unused)) ) \
            ) \
        , BOOST_MPL_PP_TOKEN_EQUAL_I \
        , 0 BOOST_PP_TUPLE_EAT(2) \
        )(a, b) \
/**/

#define BOOST_MPL_PP_TOKEN_EQUAL_I(a, b) \
    BOOST_PP_COMPL(BOOST_MPL_PP_IS_SEQ( \
        BOOST_MPL_PP_TOKEN_EQUAL_ ## a( \
            BOOST_MPL_PP_TOKEN_EQUAL_ ## b \
            )((unused)) \
        )) \
/**/

#endif // BOOST_MPL_AUX_PREPROCESSOR_TOKEN_EQUAL_HPP_INCLUDED
