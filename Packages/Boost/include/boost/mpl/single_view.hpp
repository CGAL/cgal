//-----------------------------------------------------------------------------
// boost mpl/single_view.hpp header file
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

#ifndef BOOST_MPL_SINGLE_VIEW_HPP_INCLUDED
#define BOOST_MPL_SINGLE_VIEW_HPP_INCLUDED

#include "boost/mpl/aux_/single_element_iter.hpp"
#include "boost/mpl/iterator_range.hpp"
#include "boost/mpl/aux_/void_spec.hpp"

namespace boost {
namespace mpl {

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    >
struct single_view
    : iterator_range<
          aux::single_element_iter<T,0>
        , aux::single_element_iter<T,1>
        >
{
};

BOOST_MPL_AUX_VOID_SPEC(1, single_view)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_SINGLE_VIEW_HPP_INCLUDED
