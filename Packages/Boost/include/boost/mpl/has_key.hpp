
#ifndef BOOST_MPL_HAS_KEY_HPP_INCLUDED
#define BOOST_MPL_HAS_KEY_HPP_INCLUDED

// + file: boost/mpl/has_key.hpp
// + last modified: 02/may/03

// Copyright (c) 2002-03
// David Abrahams, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.
//
// See http://www.boost.org/libs/mpl for documentation.

#include "boost/mpl/has_key_fwd.hpp"
#include "boost/mpl/aux_/has_key_impl.hpp"
#include "boost/mpl/aux_/sequence_tag.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"

namespace boost {
namespace mpl {

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(AssociativeSequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Key)
    >
struct has_key
    : has_key_impl< typename BOOST_MPL_AUX_SEQUENCE_TAG(AssociativeSequence) >
        ::template apply<AssociativeSequence,Key>
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(2,has_key,(AssociativeSequence,Key))
};

BOOST_MPL_AUX_VOID_SPEC(2, has_key)

}}

#endif // BOOST_MPL_HAS_KEY_HPP_INCLUDED
