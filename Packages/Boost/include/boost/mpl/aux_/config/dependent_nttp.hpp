//-----------------------------------------------------------------------------
// boost mpl/aux_/config/dependent_nttp.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_CONFIG_DEPENDENT_NTTP_HPP_INCLUDED
#define BOOST_MPL_AUX_CONFIG_DEPENDENT_NTTP_HPP_INCLUDED

#include "boost/config.hpp"

// GCC and EDG-based compilers incorrectly reject the following code:
//   template< typename T, T n > struct a;
//   template< typename T > struct b;
//   template< typename T, T n > struct b< a<T,n> > {};

#if    !defined(BOOST_NO_DEPENDENT_NON_TYPE_PARAMETER_IN_PARTIAL_SPECIALIZATION) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE) \
    && (   defined(__EDG__) && (__EDG_VERSION__ <= 300 || !defined(BOOST_STRICT_CONFIG)) \
        || defined(__GNUC__) && !defined(__EDG_VERSION__) && (__GNUC__ < 3 || __GNUC__ == 3 && __GNUC_MINOR__ <= 2 \
        || !defined(BOOST_STRICT_CONFIG)) \
        )

#   define BOOST_NO_DEPENDENT_NON_TYPE_PARAMETER_IN_PARTIAL_SPECIALIZATION

#endif

#endif // BOOST_MPL_AUX_CONFIG_DEPENDENT_NTTP_HPP_INCLUDED
