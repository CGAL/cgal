//-----------------------------------------------------------------------------
// boost mpl/copy_if.hpp header file
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

#ifndef BOOST_MPL_COPY_IF_HPP_INCLUDED
#define BOOST_MPL_COPY_IF_HPP_INCLUDED

#include "boost/mpl/fold.hpp"
#include "boost/mpl/aux_/copy_if_op.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/protect.hpp"
#include "boost/mpl/aux_/void_spec.hpp"

namespace boost {
namespace mpl {

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(State)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(BinaryOp)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Predicate)
    >
struct copy_if
{
 private:
    typedef typename lambda<BinaryOp>::type op_;
    typedef typename lambda<Predicate>::type pred_;

 public:
    typedef typename fold<
          Sequence
        , State
        , protect< aux::copy_if_op<op_,pred_> >
        >::type type;
};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(4, copy_if)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_COPY_IF_HPP_INCLUDED
