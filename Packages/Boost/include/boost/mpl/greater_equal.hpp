
#ifndef BOOST_MPL_GREATER_EQUAL_HPP_INCLUDED
#define BOOST_MPL_GREATER_EQUAL_HPP_INCLUDED

// Copyright (c) 2000-03 Aleksey Gurtovoy
//
// Use, modification and distribution are subject to the Boost Software 
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy 
// at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Source$
// $Date$
// $Revision$

#include "boost/mpl/bool.hpp"
#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/aux_/value_wknd.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"
#include "boost/config.hpp"

namespace boost {
namespace mpl {

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T2)
    >
struct greater_equal
{
    enum
    { 
        msvc71_wknd_ = ( BOOST_MPL_AUX_VALUE_WKND(T1)::value 
                        >= BOOST_MPL_AUX_VALUE_WKND(T2)::value )
    };

    BOOST_STATIC_CONSTANT(bool, value = msvc71_wknd_);

#if !defined(__BORLANDC__)
    typedef bool_<value> type;
#else
    typedef bool_<(
          BOOST_MPL_AUX_VALUE_WKND(T1)::value 
            >= BOOST_MPL_AUX_VALUE_WKND(T2)::value
        )> type;
#endif

    BOOST_MPL_AUX_LAMBDA_SUPPORT(2,greater_equal,(T1,T2))
};

BOOST_MPL_AUX_VOID_SPEC(2, greater_equal)

template< long N >
struct ge
{
    template< typename T > struct apply
#if !defined(__BORLANDC__)
        : greater_equal< T,integral_c<long,N> >
    {
#else
    {
        typedef typename greater_equal< T,integral_c<long,N> >::type type;
#endif
    };
};

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_GREATER_EQUAL_HPP_INCLUDED
