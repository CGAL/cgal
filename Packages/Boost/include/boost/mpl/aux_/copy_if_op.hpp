//-----------------------------------------------------------------------------
// boost mpl/aux_/copy_if_op.hpp header file
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

#ifndef BOOST_MPL_AUX_COPY_IF_OP_HPP_INCLUDED
#define BOOST_MPL_AUX_COPY_IF_OP_HPP_INCLUDED

#include "boost/mpl/identity.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/apply_if.hpp"
#include "boost/mpl/aux_/lambda_spec.hpp"

namespace boost {
namespace mpl {

namespace aux {

template<
      typename Operation
    , typename Predicate
    >
struct copy_if_op
{
    template< typename Sequence, typename T > struct apply
    {
        typedef typename apply_if<
              typename apply1<Predicate,T>::type
            , apply2<Operation,Sequence,T>
            , identity<Sequence>
            >::type type;
    };
};

} // namespace aux

BOOST_MPL_AUX_PASS_THROUGH_LAMBDA_SPEC(2,aux::copy_if_op)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_COPY_IF_OP_HPP_INCLUDED
