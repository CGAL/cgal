//-----------------------------------------------------------------------------
// boost mpl/aux_/config/ctps.hpp header file
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

#ifndef BOOST_MPL_AUX_CONFIG_CTPS_HPP_INCLUDED
#define BOOST_MPL_AUX_CONFIG_CTPS_HPP_INCLUDED

#include "boost/config.hpp"

#if    !defined(BOOST_NO_NON_TYPE_TEMPLATE_PARTIAL_SPECIALIZATION) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE) \
    && defined(__BORLANDC__) && (__BORLANDC__ <= 0x561 || !defined(BOOST_STRICT_CONFIG)) 

#   define BOOST_NO_NON_TYPE_TEMPLATE_PARTIAL_SPECIALIZATION

#endif

// BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION is defined in "boost/config.hpp"

#endif // BOOST_MPL_AUX_CONFIG_CTPS_HPP_INCLUDED
