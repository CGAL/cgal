//-----------------------------------------------------------------------------
// boost mpl/fold_backward.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2001-02
// Aleksey Gurtovoy, David Abrahams
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_FOLD_BACKWARD_HPP_INCLUDED
#define BOOST_MPL_FOLD_BACKWARD_HPP_INCLUDED

#include "boost/mpl/begin_end.hpp"
#include "boost/mpl/O1_size.hpp"
#include "boost/mpl/arg.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/aux_/fold_backward_impl.hpp"
#include "boost/mpl/aux_/void_spec.hpp"

namespace boost {
namespace mpl {

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(State)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(BackwardOp)
    , typename ForwardOp = arg<1>
    >
struct fold_backward
{
    typedef typename aux::fold_backward_impl<
          ::boost::mpl::O1_size<Sequence>::value
        , typename begin<Sequence>::type
        , typename end<Sequence>::type
        , State
        , typename lambda<BackwardOp>::type
        , typename lambda<ForwardOp>::type
        >::state type;

    BOOST_MPL_AUX_LAMBDA_SUPPORT(3,fold_backward,(Sequence,State,BackwardOp))
};

BOOST_MPL_AUX_VOID_SPEC(3, fold_backward)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_FOLD_BACKWARD_HPP_INCLUDED
