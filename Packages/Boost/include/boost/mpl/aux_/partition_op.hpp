//-----------------------------------------------------------------------------
// boost mpl/aux_/partition_op.hpp header file
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

#ifndef BOOST_MPL_AUX_PARTITION_OP_HPP_INCLUDED
#define BOOST_MPL_AUX_PARTITION_OP_HPP_INCLUDED

#include "boost/mpl/apply.hpp"
#include "boost/mpl/apply_if.hpp"
#include "boost/mpl/if.hpp"
#include "boost/mpl/pair.hpp"
#include "boost/mpl/push_front.hpp"
#include "boost/mpl/aux_/lambda_spec.hpp"

namespace boost {
namespace mpl {
namespace aux {

template <typename Predicate>
struct partition_op
{
    template <typename State, typename Iter>
    struct apply
    {
    private:
        typedef typename State::first first_;
        typedef typename State::second second_;
        typedef typename Iter::type t_;
        typedef typename apply1< Predicate,t_ >::type pred_;

        typedef typename apply_if<
              pred_
            , push_front< first_, t_ >
            , push_front< second_, t_ >
            >::type result_;

    public:
        typedef typename if_<
              pred_
            , pair< result_,second_ >
            , pair< first_,result_ >
            >::type type;
    };
};

} // namespace aux

BOOST_MPL_AUX_PASS_THROUGH_LAMBDA_SPEC(1,aux::partition_op)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_PARTITION_OP_HPP_INCLUDED
