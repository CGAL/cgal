//-----------------------------------------------------------------------------
// boost mpl/aux_/config/bind.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002
// David Abrahams, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_CONFIG_BIND_HPP_INCLUDED
#define BOOST_MPL_AUX_CONFIG_BIND_HPP_INCLUDED

#include "boost/config.hpp"

#if    !defined(BOOST_MPL_NO_BIND_TEMPLATE) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE) \
    && (   defined(BOOST_MSVC) && BOOST_MSVC <= 1300 \
        || defined(__BORLANDC__) && (__BORLANDC__ <= 0x561 || !defined(BOOST_STRICT_CONFIG)) \
        )

#   define BOOST_MPL_NO_BIND_TEMPLATE

#endif 

#endif // BOOST_MPL_AUX_CONFIG_BIND_HPP_INCLUDED
