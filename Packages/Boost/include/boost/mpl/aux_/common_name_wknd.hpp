//-----------------------------------------------------------------------------
// boost mpl/aux_/common_name_wknd.hpp header file
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

#ifndef BOOST_MPL_AUX_COMMON_NAME_WKND_HPP_INCLUDED
#define BOOST_MPL_AUX_COMMON_NAME_WKND_HPP_INCLUDED

#if defined(__BORLANDC__) && __BORLANDC__ < 0x561
// agurt 12/nov/02: to suppress the bogus "Cannot have both a template class and 
// function named 'xxx'" diagnostic
#   define BOOST_MPL_AUX_COMMON_NAME_WKND(name) \
namespace name_##wknd { \
template< typename > void name(); \
} \
/**/

#else

#   define BOOST_MPL_AUX_COMMON_NAME_WKND(name) /**/

#endif // __BORLANDC__

#endif // BOOST_MPL_AUX_COMMON_NAME_WKND_HPP_INCLUDED
