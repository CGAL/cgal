//-----------------------------------------------------------------------------
// boost mpl/vector/aux_/item.hpp header file
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

#ifndef BOOST_MPL_VECTOR_AUX_ITEM_HPP_INCLUDED
#define BOOST_MPL_VECTOR_AUX_ITEM_HPP_INCLUDED

#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/aux_/config/vector.hpp"
#include "boost/mpl/aux_/config/ctps.hpp"

namespace boost {
namespace mpl {

#if defined(BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL)

template< typename V, long N >
struct vector_item
{
    typedef __typeof__(V::item(integral_c<long,V::size::value - N - 1>())) wrapped_type_;
    typedef typename wrapped_type_::type type;
};

#else

// agurt, 29/sep/02: Borland 5.5.1 does a poor job on partial 
// specialization here
#   if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) \
    && !defined(BOOST_NO_NON_TYPE_TEMPLATE_PARTIAL_SPECIALIZATION)

// forward declaration
template< typename T, long N >
struct vector_item;

#   else

namespace aux {

template< long N > struct vector_item_impl
{
    template< typename V > struct result_;
};

} // namespace aux

template< typename T, long N >
struct vector_item
    : aux::vector_item_impl<N>::template result_<T>
{
};

#   endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#endif // BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_VECTOR_AUX_ITEM_HPP_INCLUDED
