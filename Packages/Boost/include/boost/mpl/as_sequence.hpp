//-----------------------------------------------------------------------------
// boost mpl/as_sequence.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AS_SEQUENCE_HPP_INCLUDED
#define BOOST_MPL_AS_SEQUENCE_HPP_INCLUDED

#include "boost/mpl/is_sequence.hpp"
#include "boost/mpl/single_view.hpp"
#include "boost/mpl/if.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"
#include "boost/mpl/aux_/config/eti.hpp"

namespace boost { namespace mpl {

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    >
struct as_sequence
    : if_< is_sequence<T>, T, single_view<T> >
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1,as_sequence,(T))
};

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
template<> struct as_sequence<int>
{
    typedef single_view<int> type;
};
#endif

BOOST_MPL_AUX_VOID_SPEC(1, as_sequence)

}} // namespace boost::mpl

#endif // BOOST_MPL_AS_SEQUENCE_HPP_INCLUDED
