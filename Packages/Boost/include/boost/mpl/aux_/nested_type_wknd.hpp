//-----------------------------------------------------------------------------
// boost mpl/aux_/nested_type_wknd.hpp header file
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

#ifndef BOOST_MPL_AUX_NESTED_TYPE_WKND_HPP_INCLUDED
#define BOOST_MPL_AUX_NESTED_TYPE_WKND_HPP_INCLUDED

#include "boost/config.hpp"

#if defined(__GNUC__) && (__GNUC__ < 3 || __GNUC__ == 3 && __GNUC_MINOR__ <= 2 \
    || !defined(BOOST_STRICT_CONFIG)) \
 || defined(__BORLANDC__) && (__BORLANDC__ <= 0x561 || !defined(BOOST_STRICT_CONFIG)) \
 || defined(__SUNPRO_CC)

namespace boost { namespace mpl { namespace aux {

template< typename T >
struct nested_type_wknd
    : T::type
{
};

}}} // namespace boost::mpl::aux

#   define BOOST_MPL_AUX_NESTED_TYPE_WKND(T) ::boost::mpl::aux::nested_type_wknd<T>

#else

#   define BOOST_MPL_AUX_NESTED_TYPE_WKND(T) T::type

#endif // __GNUC__

#endif // BOOST_MPL_AUX_NESTED_TYPE_WKND_HPP_INCLUDED
