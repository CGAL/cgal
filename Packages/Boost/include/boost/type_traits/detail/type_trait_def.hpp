//-----------------------------------------------------------------------------
// boost/type_traits/detail/type_trait_def.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002
// Aleksey Gurtovoy
//
// Use, modification and distribution are subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

// no include guards, the header is intended for multiple inclusion!

#include "boost/type_traits/detail/template_arity_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"

#define BOOST_TT_AUX_TYPE_TRAIT_DEF1(trait,T,result) \
template< typename T > struct trait \
{ \
    typedef result type; \
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1,trait,(T)) \
}; \
\
BOOST_TT_AUX_TEMPLATE_ARITY_SPEC(1,trait) \
/**/

#define BOOST_TT_AUX_TYPE_TRAIT_SPEC1(trait,spec,result) \
template<> struct trait<spec> \
{ \
    typedef result type; \
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(1,trait,(spec)) \
}; \
/**/

#define BOOST_TT_AUX_TYPE_TRAIT_IMPL_SPEC1(trait,spec,result) \
template<> struct trait##_impl<spec> \
{ \
    typedef result type; \
}; \
/**/

#define BOOST_TT_AUX_TYPE_TRAIT_PARTIAL_SPEC1_1(param,trait,spec,result) \
template< param > struct trait<spec> \
{ \
    typedef result type; \
}; \
/**/

#define BOOST_TT_AUX_TYPE_TRAIT_PARTIAL_SPEC1_2(param1,param2,trait,spec,result) \
template< param1, param2 > struct trait<spec> \
{ \
    typedef result; \
}; \
/**/

#define BOOST_TT_AUX_TYPE_TRAIT_IMPL_PARTIAL_SPEC1_1(param,trait,spec,result) \
template< param > struct trait##_impl<spec> \
{ \
    typedef result type; \
}; \
/**/
