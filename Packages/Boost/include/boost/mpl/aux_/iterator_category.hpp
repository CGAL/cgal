//-----------------------------------------------------------------------------
// boost mpl/aux_/iterator_category.hpp header file
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

#ifndef BOOST_MPL_AUX_ITERATOR_CATEGORY_HPP_INCLUDED
#define BOOST_MPL_AUX_ITERATOR_CATEGORY_HPP_INCLUDED

#include "boost/mpl/aux_/config/internals.hpp"

#if defined(BOOST_MPL_INTERNALS_USE_ITERATOR_CATEGORY)

#   include "boost/mpl/iterator_category.hpp"
#   define BOOST_MPL_AUX_ITERATOR_CATEGORY(iter) iterator_category<iter>::type

#else

#   define BOOST_MPL_AUX_ITERATOR_CATEGORY(iter) iter::category

#endif

#endif // BOOST_MPL_AUX_ITERATOR_CATEGORY_HPP_INCLUDED
