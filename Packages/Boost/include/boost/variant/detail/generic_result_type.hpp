//-----------------------------------------------------------------------------
// boost variant/detail/generic_result_type.hpp header file
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

#ifndef BOOST_VARIANT_DETAIL_GENERIC_RESULT_TYPE_HPP
#define BOOST_VARIANT_DETAIL_GENERIC_RESULT_TYPE_HPP

#include "boost/config.hpp"

//////////////////////////////////////////////////////////////////////////
// (workaround) macro BOOST_VARIANT_AUX_GENERIC_RESULT_TYPE
//
// On compilers with BOOST_NO_VOID_RETURNS, this macro provides a route
// to a single syntax for dealing with template functions that may (but
// not necessarily) return nothing (i.e. void).
//
// BOOST_VARIANT_AUX_RETURN_VOID provided for compilers w/ (erroneous?)
// warnings about non-void functions not returning a value.
//

#if !defined(BOOST_NO_VOID_RETURNS)

#define BOOST_VARIANT_AUX_GENERIC_RESULT_TYPE(T) \
    T   \
    /**/

#define BOOST_VARIANT_AUX_RETURN_VOID  \
    /**/

#define BOOST_VARIANT_AUX_RETURN_VOID_TYPE \
    void    \
    /**/

#else // defined(BOOST_NO_VOID_RETURNS)

namespace boost {
namespace detail { namespace variant {

struct fake_return_void
{
    fake_return_void()
    {
    }

    template <typename T>
    fake_return_void(const T&)
    {
    }
};

template <typename T>
struct no_void_returns_helper
{
    typedef T type;
};

template <>
struct no_void_returns_helper<void>
{
    typedef fake_return_void type;
};

}} // namespace detail::variant
} // namespace boost

#define BOOST_VARIANT_AUX_GENERIC_RESULT_TYPE(T) \
    BOOST_DEDUCED_TYPENAME                                           \
        ::boost::detail::variant::no_void_returns_helper< T >::type  \
    /**/

#define BOOST_VARIANT_AUX_RETURN_VOID  \
    return ::boost::detail::variant::fake_return_void()     \
    /**/

#define BOOST_VARIANT_AUX_RETURN_VOID_TYPE  \
    ::boost::detail::variant::fake_return_void

#endif // BOOST_NO_VOID_RETURNS workaround

#endif // BOOST_VARIANT_DETAIL_GENERIC_RESULT_TYPE_HPP
