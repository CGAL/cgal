//-----------------------------------------------------------------------------
// boost mpl/aux_/O1_size_impl.hpp header file
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

#ifndef BOOST_MPL_O1_SIZE_IMPL_HPP_INCLUDED
#define BOOST_MPL_O1_SIZE_IMPL_HPP_INCLUDED

#include "boost/mpl/O1_size_fwd.hpp"
#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/if.hpp"
#include "boost/mpl/aux_/has_size.hpp"

namespace boost {
namespace mpl {

// default implementation - returns |Sequence::size| if sequence has a |size|
// member, and -1 otherwise; conrete sequences might override it by specializing 
// either the |O1_size_traits| or the primary |O1_size| template

#if 0//!defined(BOOST_MSVC) || BOOST_MSVC > 1300
namespace aux {

template< typename Sequence >
struct O1_size_impl
    : Sequence::size
{
};

} // namespace aux

template< typename Tag >
struct O1_size_traits
{
    template< typename Sequence > struct algorithm
        : if_c<
              ::boost::mpl::aux::has_size<Sequence>::value
            , aux::O1_size_impl<Sequence>
            , integral_c<long,-1>
            >::type
    {
    };
};
#else

template< typename Tag >
struct O1_size_traits
{
    template< typename Sequence > struct algorithm
        : integral_c<long,-1>
        {
        };
};

#endif // BOOST_MSVC > 1300


} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_O1_SIZE_IMPL_HPP_INCLUDED
