
// NO INCLUDE GUARDS, THE HEADER IS INTENDED FOR MULTIPLE INCLUSION

// Copyright Aleksey Gurtovoy 2002-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)

// $Source$
// $Date$
// $Revision$

#include <boost/type_traits/detail/template_arity_spec.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/mpl/aux_/lambda_support.hpp>

#include <cstddef>

#if !defined(BOOST_MSVC) || BOOST_MSVC > 1200
#   define BOOST_TT_AUX_SIZE_T_TRAIT_VALUE_DECL(C) /**/
#else
#   define BOOST_TT_AUX_SIZE_T_TRAIT_VALUE_DECL(C) \
    typedef mpl::size_t< C > base_; \
    using base_::value; \
    /**/
#endif


#define BOOST_TT_AUX_SIZE_T_TRAIT_DEF1(trait,T,C) \
template< typename T > struct trait \
    : mpl::size_t< C > \
{ \
    BOOST_TT_AUX_SIZE_T_TRAIT_VALUE_DECL(C) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1,trait,(T)) \
}; \
\
BOOST_TT_AUX_TEMPLATE_ARITY_SPEC(1,trait) \
/**/

#define BOOST_TT_AUX_SIZE_T_TRAIT_SPEC1(trait,spec,C) \
template<> struct trait<spec> \
    : mpl::size_t< C > \
{ \
    BOOST_TT_AUX_SIZE_T_TRAIT_VALUE_DECL(C) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(1,trait,(spec)) \
}; \
/**/

#define BOOST_TT_AUX_SIZE_T_TRAIT_PARTIAL_SPEC1_1(param,trait,spec,C) \
template< param > struct trait<spec> \
    : mpl::size_t< C > \
{ \
}; \
/**/
