//-----------------------------------------------------------------------------
// boost mpl/empty_base.hpp header file
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

#ifndef BOOST_MPL_EMPTY_BASE_HPP_INCLUDED
#define BOOST_MPL_EMPTY_BASE_HPP_INCLUDED

#include "boost/type_traits/is_empty.hpp"
#include "boost/mpl/bool.hpp"
#include "boost/config.hpp"

// should be always the last #include directive
#include "boost/type_traits/detail/bool_trait_def.hpp"

namespace boost {
namespace mpl {

// empty base class, guaranteed to have no members; inheritance from
// 'empty_base' through the 'inherit' metafunction is a no-op - see 
// "mpl/inherit.hpp" header for the details
struct empty_base {};

template< typename T >
struct is_empty_base
    : false_
{
#if defined(BOOST_MSVC) && BOOST_MSVC < 1300
    using false_::value;
#endif
};

template<>
struct is_empty_base<empty_base>
    : true_
{
#if defined(BOOST_MSVC) && BOOST_MSVC < 1300
    using true_::value;
#endif
};

} // namespace mpl

BOOST_TT_AUX_BOOL_TRAIT_SPEC1(is_empty, mpl::empty_base, true)

} // namespace boost

#include "boost/type_traits/detail/bool_trait_undef.hpp"

#endif // BOOST_MPL_EMPTY_BASE_HPP_INCLUDED
