//-----------------------------------------------------------------------------
// boost mpl/aux_/config/vector.hpp header file
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

#ifndef BOOST_MPL_AUX_CONFIG_VECTOR_HPP_INCLUDED
#define BOOST_MPL_AUX_CONFIG_VECTOR_HPP_INCLUDED

// agurt, 10/jul/02: full-fledged __typeof is needed to permit the optimal 
// vector implementation

#if    !defined(BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE) \
    && defined(__MWERKS__) && __MWERKS__ >= 0x3001
    
#   define BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL

#endif

#endif // BOOST_MPL_AUX_CONFIG_VECTOR_HPP_INCLUDED
