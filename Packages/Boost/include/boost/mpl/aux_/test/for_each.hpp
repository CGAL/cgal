
#ifndef BOOST_MPL_AUX_TEST_FOR_EACH_HPP_INCLUDED
#define BOOST_MPL_AUX_TEST_FOR_EACH_HPP_INCLUDED

// + file: boost/mpl/aux_/test/for_each.hpp
// + last modified: 04/may/03

// Copyright (c) 2002-03
// Aleksey Gurtovoy
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

#include "boost/preprocessor/seq/for_each_i.hpp"
#include "boost/preprocessor/seq/transform.hpp" 
#include "boost/preprocessor/seq/enum.hpp" 
#include "boost/preprocessor/tuple/elem.hpp"
#include "boost/preprocessor/control/expr_if.hpp"
#include "boost/preprocessor/cat.hpp"

#define CTT_AUX_IS_PLACEHOLDER_1 (1)

#define CTT_AUX_IS_PLACEHOLDER(x) \
    BOOST_MPL_PP_IS_SEQ( BOOST_PP_CAT(CTT_AUX_IS_PLACEHOLDER, x) ) \
/**/

#define CTT_AUX_PARAM_OR_T( unused, T, param ) \
    BOOST_PP_TUPLE_ELEM( 2, CTT_AUX_IS_PLACEHOLDER(param), (param, T) ) \
/**/

#define CTT_AUX_REPLACE_PLACEHOLDER( params_seq, T ) \
    BOOST_PP_SEQ_ENUM( \
          BOOST_PP_SEQ_TRANSFORM( CTT_AUX_PARAM_OR_T, T, params_seq ) \
        ) \
/**/

#define CTT_AUX_INVOKE_TEST( unused, test_params_pair, i, T ) \
    BOOST_PP_EXPR_IF(i, ;) \
    BOOST_PP_TUPLE_ELEM(2, 0, test_params_pair)< \
        CTT_AUX_REPLACE_PLACEHOLDER( \
              BOOST_PP_TUPLE_ELEM(2, 1, test_params_pair) \
            , T \
            ) \
        >() \
/**/

#define CTT_for_each( seq, test, params_seq ) \
    BOOST_PP_SEQ_FOR_EACH_I( CTT_AUX_INVOKE_TEST, (test, params_seq), seq ) \
/**/

#endif // BOOST_MPL_AUX_TEST_FOR_EACH_HPP_INCLUDED
