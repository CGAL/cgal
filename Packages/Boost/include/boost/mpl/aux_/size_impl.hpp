//-----------------------------------------------------------------------------
// boost mpl/aux_/size_impl.hpp header file
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

#ifndef BOOST_MPL_AUX_SIZE_IMPL_HPP_INCLUDED
#define BOOST_MPL_AUX_SIZE_IMPL_HPP_INCLUDED

#include "boost/mpl/size_fwd.hpp"
#include "boost/mpl/begin_end.hpp"
#include "boost/mpl/distance.hpp"
#include "boost/mpl/aux_/traits_lambda_spec.hpp"
#include "boost/config.hpp"

namespace boost {
namespace mpl {

// default implementation; conrete sequences might override it by 
// specializing either the 'size_traits' or the primary 'size' template

template< typename Tag >
struct size_traits
{
    template< typename Sequence > struct algorithm
#if !defined(__BORLANDC__) && (__BORLANDC__ <= 0x561 || !defined(BOOST_STRICT_CONFIG))
        : distance<
              typename begin<Sequence>::type
            , typename end<Sequence>::type
            >
    {
#else
    {
        typedef typename distance<
              typename begin<Sequence>::type
            , typename end<Sequence>::type
            >::type type;
#endif
    };
};

BOOST_MPL_ALGORITM_TRAITS_LAMBDA_SPEC(1,size_traits)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_SIZE_IMPL_HPP_INCLUDED
