//-----------------------------------------------------------------------------
// boost mpl/aux_/sequence_tag.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2001-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_SEQUENCE_TAG_HPP_INCLUDED
#define BOOST_MPL_AUX_SEQUENCE_TAG_HPP_INCLUDED

//#include "boost/mpl/aux_/config/internal.hpp"
#include "boost/config.hpp"

#if !defined(BOOST_MPL_INTERNAL_USE_SEQUENCE_TAG)
#   define BOOST_MPL_INTERNAL_USE_SEQUENCE_TAG
#endif

#if defined(BOOST_MPL_INTERNAL_USE_SEQUENCE_TAG) \
 || defined(BOOST_MSVC) && BOOST_MSVC < 1300

#   include "boost/mpl/sequence_tag.hpp"
#   define BOOST_MPL_AUX_SEQUENCE_TAG(seq) sequence_tag<seq>::type

#else

#   define BOOST_MPL_AUX_SEQUENCE_TAG(seq) seq::tag

#endif

#endif // BOOST_MPL_AUX_SEQUENCE_TAG_HPP_INCLUDED
