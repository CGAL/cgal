//-----------------------------------------------------------------------------
// boost mpl/push_back.hpp header file
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

#ifndef BOOST_MPL_PUSH_BACK_HPP_INCLUDED
#define BOOST_MPL_PUSH_BACK_HPP_INCLUDED

#include "boost/mpl/push_back_fwd.hpp"
#include "boost/mpl/aux_/push_back_impl.hpp"
#include "boost/mpl/aux_/sequence_tag.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"

namespace boost {
namespace mpl {

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    >
struct push_back
    : push_back_traits< typename BOOST_MPL_AUX_SEQUENCE_TAG(Sequence) >
        ::template algorithm< Sequence,T >
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(2,push_back,(Sequence,T))
};

BOOST_MPL_AUX_VOID_SPEC(2, push_back)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_PUSH_BACK_HPP_INCLUDED
