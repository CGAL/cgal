//-----------------------------------------------------------------------------
// boost mpl/aux_/algorithm_namespace.hpp header file
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

#ifndef BOOST_MPL_AUX_ALGORITHM_NAMESPACE_HPP_INCLUDED
#define BOOST_MPL_AUX_ALGORITHM_NAMESPACE_HPP_INCLUDED

#if defined(__GNUC__) && __GNUC__ >= 3

#   define BOOST_MPL_AUX_AGLORITHM_NAMESPACE_PREFIX algo_::
#   define BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN namespace algo_ {
#   define BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END } using namespace algo_;

#else

#   define BOOST_MPL_AUX_AGLORITHM_NAMESPACE_PREFIX /**/
#   define BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN /**/
#   define BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END /**/

#endif

#endif // BOOST_MPL_AUX_ALGORITHM_NAMESPACE_HPP_INCLUDED
