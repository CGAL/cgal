
#ifndef BOOST_MPL_AUX_TEST_ASSERT_HPP_INCLUDED
#define BOOST_MPL_AUX_TEST_ASSERT_HPP_INCLUDED

// + file: boost/mpl/aux_/test/assert.hpp
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

#include "boost/mpl/aux_/config/ctps.hpp"
#include "boost/mpl/aux_/config/msvc.hpp"
#include "boost/mpl/aux_/config/workaround.hpp"

#include "boost/static_assert.hpp"

#include "boost/preprocessor/tuple/rem.hpp"
#include "boost/preprocessor/control/expr_if.hpp"
#include "boost/preprocessor/logical/not.hpp"

#define CTT_assert( expr ) BOOST_STATIC_ASSERT( expr )
#define CTT_assert_equal(arity, tuple) assert_equal< BOOST_PP_TUPLE_REM(arity)tuple >::type()
#define CTT_assert_not_equal(arity, tuple) assert_not_equal< BOOST_PP_TUPLE_REM(arity)tuple >()
#define CTT_assert_same(arity, tuple) assert_same< BOOST_PP_TUPLE_REM(arity)tuple >()
#define CTT_assert_not_same(arity, tuple) assert_not_same< BOOST_PP_TUPLE_REM(arity)tuple >()

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
#   define AUX_ASSERT_DEF( param_type, suffix, test_equality ) \
template< param_type x, param_type y > struct assert_##suffix \
BOOST_PP_EXPR_IF( BOOST_PP_NOT(test_equality), {} ); \
template< param_type x > struct assert_##suffix<x,x> \
BOOST_PP_EXPR_IF( test_equality, {} ); \
/**/
#elif BOOST_WORKAROUND(BOOST_MSVC, <= 1300)
#   define AUX_ASSERT_DEF( param_type, suffix, test_equality ) \
template< param_type x > struct assert_##suffix##_impl \
{ \
    template< param_type y > struct result_; \
    template<> struct result_<x> \
    { BOOST_PP_EXPR_IF( BOOST_PP_NOT(test_equality), private: virtual ~result_<x>() = 0; ) }; \
    template< param_type y > struct result_ \
    { BOOST_PP_EXPR_IF( test_equality, private: virtual ~result_<y>() = 0;  ) }; \
}; \
\
template< param_type x, param_type y > struct assert_##suffix \
    : assert_##suffix##_impl<x>::template result_<y> {}; \
/**/
#else
#   define AUX_ASSERT_DEF( param_type, suffix, test_equality ) \
template< param_type x > struct assert_##suffix##_impl \
{ \
    template< param_type y > struct result_ \
    BOOST_PP_EXPR_IF( BOOST_PP_NOT(test_equality), {} ); \
}; \
template< param_type x > struct assert_##suffix##_impl<x>::result_<x> \
BOOST_PP_EXPR_IF( test_equality, {} ); \
\
template< param_type x, param_type y > struct assert_##suffix \
    : assert_##suffix##_impl<x>::template result_<y> {}; \
/**/
#endif

namespace {
AUX_ASSERT_DEF(long, equal, 1)
AUX_ASSERT_DEF(long, not_equal, 0)
AUX_ASSERT_DEF(typename, same, 1)
AUX_ASSERT_DEF(typename, not_same, 0)
}

#undef AUX_ASSERT_DEF

#endif // BOOST_MPL_AUX_TEST_ASSERT_HPP_INCLUDED
