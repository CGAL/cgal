//-----------------------------------------------------------------------------
// boost mpl/aux_/value_wknd.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_VALUE_WKND_HPP_INCLUDED
#define BOOST_MPL_AUX_VALUE_WKND_HPP_INCLUDED

#include "boost/mpl/aux_/config/eti.hpp"

#if defined(__BORLANDC__) && (__BORLANDC__ <= 0x561 || !defined(BOOST_STRICT_CONFIG)) \
 || defined(BOOST_MPL_MSVC_60_ETI_BUG)
 
#   include "boost/mpl/int.hpp"

namespace boost { namespace mpl { namespace aux {

template< typename C_ >
struct value_wknd
    : C_
{
};

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
template<>
struct value_wknd<int>
    : int_<1>
{
};
#endif

}}} // namespace boost::mpl::aux

#   if !defined(BOOST_MPL_MSVC_60_ETI_BUG)
#       define BOOST_MPL_AUX_VALUE_WKND(C) ::boost::mpl::aux::value_wknd< C >
#       define BOOST_MPL_AUX_MSVC_VALUE_WKND(C) BOOST_MPL_AUX_VALUE_WKND(C)
#   else
#       define BOOST_MPL_AUX_VALUE_WKND(C) C
#       define BOOST_MPL_AUX_MSVC_VALUE_WKND(C) ::boost::mpl::aux::value_wknd< C >
#   endif

#else

#   define BOOST_MPL_AUX_VALUE_WKND(C) C
#   define BOOST_MPL_AUX_MSVC_VALUE_WKND(C) C

#endif // __BORLANDC__ || BOOST_MPL_MSVC_60_ETI_BUG

#endif // BOOST_MPL_AUX_VALUE_WKND_HPP_INCLUDED
