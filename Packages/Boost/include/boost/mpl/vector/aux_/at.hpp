//-----------------------------------------------------------------------------
// boost mpl/vector/aux_/at.hpp header file
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

#ifndef BOOST_MPL_VECTOR_AUX_AT_HPP_INCLUDED
#define BOOST_MPL_VECTOR_AUX_AT_HPP_INCLUDED

#include "boost/mpl/at_fwd.hpp"
#include "boost/mpl/aux_/value_wknd.hpp"
#include "boost/mpl/vector/aux_/item.hpp"
#include "boost/mpl/vector/aux_/tag.hpp"
#include "boost/mpl/aux_/config/vector.hpp"
#include "boost/mpl/aux_/config/ctps.hpp"

namespace boost {
namespace mpl {

#if defined(BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL)

template<>
struct at_traits< aux::vector_tag >
{
    template< typename Vector, typename N > struct algorithm
        : vector_item<
              Vector
            , BOOST_MPL_AUX_VALUE_WKND(N)::value
            >
    {
    };
};

#else

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) \
 && !defined(BOOST_NO_NON_TYPE_TEMPLATE_PARTIAL_SPECIALIZATION)

template< long S >
struct at_traits< aux::vector_tag<S> >
{
    template< typename Vector, typename N > struct algorithm
#if !defined(__BORLANDC__)
        : vector_item<
              Vector
            , BOOST_MPL_AUX_VALUE_WKND(N)::value
            >
    {
#else
    {
        typedef typename vector_item<
              Vector
            , BOOST_MPL_AUX_VALUE_WKND(N)::value
            >::type type;
#endif
    };
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#endif // BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_VECTOR_AUX_AT_HPP_INCLUDED
