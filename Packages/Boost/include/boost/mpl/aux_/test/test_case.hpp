
#ifndef BOOST_MPL_AUX_TEST_TEST_CASE_HPP_INCLUDED
#define BOOST_MPL_AUX_TEST_TEST_CASE_HPP_INCLUDED

// + file: boost/mpl/aux_/test/test_case.hpp
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

#include "boost/preprocessor/comma_if.hpp"
#include "boost/preprocessor/seq/for_each_i.hpp"
#include "boost/preprocessor/seq/enum.hpp" 

#define CTT_AUX_PARAM_DEF(unused, prefix, i, elem) \
    BOOST_PP_COMMA_IF(i) prefix elem \
/**/

#define CTT_test_case( name, params_seq ) \
template< \
      BOOST_PP_SEQ_FOR_EACH_I( CTT_AUX_PARAM_DEF, typename, params_seq ) \
    > \
void name() \
/**/

#define CTT_test( test, params_seq ) \
    test< BOOST_PP_SEQ_ENUM(params_seq) >() \
/**/

#endif // BOOST_MPL_AUX_TEST_TEST_CASE_HPP_INCLUDED
