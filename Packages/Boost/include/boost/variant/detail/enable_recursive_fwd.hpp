//-----------------------------------------------------------------------------
// boost variant/detail/enable_recursive_fwd.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2003
// Eric Friedman
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_VARIANT_DETAIL_ENABLE_RECURSIVE_FWD_HPP
#define BOOST_VARIANT_DETAIL_ENABLE_RECURSIVE_FWD_HPP

#include "boost/mpl/aux_/config/ctps.hpp"

#include "boost/mpl/bool_fwd.hpp"

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
#   include "boost/mpl/bool.hpp"
#else
#   include "boost/type_traits/is_base_and_derived.hpp"
#endif

namespace boost {
namespace detail { namespace variant {

///////////////////////////////////////////////////////////////////////////////
// (detail) tag recursive_flag
//
// Signifies that the variant should perform recursive substituion.
//

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

template <typename T>
struct recursive_flag
{
    typedef T type;
};

#else // defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

struct recursive_flag_tag
{
};

template <typename T>
struct recursive_flag
    : recursive_flag_tag
{
    typedef T type;
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION workaround

///////////////////////////////////////////////////////////////////////////////
// (detail) metafunction is_recursive_flag
//
// Signifies that the variant should perform recursive substituion.
//

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

template <typename T>
struct is_recursive_flag
    : mpl::false_
{
};

template <typename T>
struct is_recursive_flag< recursive_flag<T> >
    : mpl::true_
{
};

#else // defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

template <typename T>
struct is_recursive_flag
    : is_base_and_derived< recursive_flag_tag,T >
{
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION workaround

///////////////////////////////////////////////////////////////////////////////
// (detail) metafunction enable_recursive
//
// Attempts recursive_variant_ tag substitution, wrapping with
// boost::recursive_wrapper if substituion occurs w/ non-indirect result
// (i.e., not a reference or pointer) *and* NoWrapper is false_.
//
template <
      typename T
    , typename RecursiveVariant
    , typename NoWrapper = mpl::false_
    >
struct enable_recursive;

///////////////////////////////////////////////////////////////////////////////
// (detail) metafunction class quoted_enable_recursive
//
// Same behavior as enable_recursive metafunction (see above).
//
template <
      typename RecursiveVariant
    , typename NoWrapper = mpl::false_
    >
struct quoted_enable_recursive;

}} // namespace detail::variant
} // namespace boost

#endif // BOOST_VARIANT_DETAIL_ENABLE_RECURSIVE_FWD_HPP
