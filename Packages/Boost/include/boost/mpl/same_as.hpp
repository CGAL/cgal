//-----------------------------------------------------------------------------
// boost mpl/same_as.hpp header file
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

#ifndef BOOST_MPL_SAME_AS_HPP_INCLUDED
#define BOOST_MPL_SAME_AS_HPP_INCLUDED

#include "boost/mpl/not.hpp"
#include "boost/mpl/aux_/lambda_spec.hpp"
#include "boost/type_traits/is_same.hpp"

namespace boost {
namespace mpl {

template< typename T1 >
struct same_as
{
    template< typename T2 > struct apply
#if !defined(__BORLANDC__) || (__BORLANDC__ > 0x551 && defined(BOOST_STRICT_CONFIG))
        : is_same<T1,T2>
    {
#else
    {
        typedef typename is_same<T1,T2>::type type;
#endif
    };
};

template< typename T1 >
struct not_same_as
{
    template< typename T2 > struct apply
#if !defined(__BORLANDC__) || (__BORLANDC__ > 0x51 && defined(BOOST_STRICT_CONFIG))
        : not_< is_same<T1,T2> >
    {
#else
    {
        typedef typename not_< is_same<T1,T2> >::type type;
#endif
    };
};

BOOST_MPL_AUX_PASS_THROUGH_LAMBDA_SPEC(1,same_as)
BOOST_MPL_AUX_PASS_THROUGH_LAMBDA_SPEC(1,not_same_as)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_SAME_AS_HPP_INCLUDED
