
#ifndef BOOST_MPL_SHIFT_RIGHT_HPP_INCLUDED
#define BOOST_MPL_SHIFT_RIGHT_HPP_INCLUDED

// + file: boost/mpl/shift_right.hpp
// + last modified: 25/feb/03

// Copyright (c) 2000-03
// Aleksey Gurtovoy, Jaap Suter
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

#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/aux_/typeof.hpp"
#include "boost/mpl/aux_/value_wknd.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"
#include "boost/config.hpp"

namespace boost {
namespace mpl {

template<
      typename T, T N, typename ShiftT, ShiftT Shift
    >
struct shift_right_c
{
    BOOST_STATIC_CONSTANT(T, value = (N >> Shift));
#if !defined(__BORLANDC__)
    typedef integral_c<T,value> type;
#else
    typedef integral_c<T,(N >> Shift)> type;
#endif
};

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Shift)
    >
struct shift_right
    : shift_right_c<
          BOOST_MPL_AUX_TYPEOF(T, T::value)
        , BOOST_MPL_AUX_MSVC_VALUE_WKND(T)::value
        , BOOST_MPL_AUX_TYPEOF(Shift, Shift::value)
        , BOOST_MPL_AUX_MSVC_VALUE_WKND(Shift)::value
        >
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(2, shift_right, (T,Shift))
};

BOOST_MPL_AUX_VOID_SPEC(2,shift_right)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_SHIFT_RIGHT_HPP_INCLUDED
