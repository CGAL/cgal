//-----------------------------------------------------------------------------
// boost mpl/sort_fwd.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002-2003
// Eric Friedman
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_SORT_FWD_HPP_INCLUDED
#define BOOST_MPL_SORT_FWD_HPP_INCLUDED

#include "boost/mpl/less.hpp"
#include "boost/mpl/placeholders.hpp"
#include "boost/mpl/aux_/void_spec.hpp"

namespace boost {
namespace mpl {

template< typename Tag > struct sort_traits;

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template <
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename Predicate = less<_,_>
    >
struct sort;

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(1, sort)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_SORT_FWD_HPP_INCLUDED
