//-----------------------------------------------------------------------------
// boost mpl/iter_fold_if.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2003
// Eric Friedman, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_ITER_FOLD_IF_HPP_INCLUDED
#define BOOST_MPL_ITER_FOLD_IF_HPP_INCLUDED

#include "boost/mpl/aux_/config/workaround.hpp"
#include "boost/mpl/aux_/iter_fold_if_impl.hpp"
#include "boost/mpl/and.hpp"
#include "boost/mpl/always.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/apply_if.hpp"
#include "boost/mpl/begin_end.hpp"
#include "boost/mpl/bool.hpp"
#include "boost/mpl/if.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/not.hpp"
#include "boost/mpl/pair.hpp"
#include "boost/mpl/void.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"
#include "boost/type_traits/is_same.hpp"

namespace boost {
namespace mpl {

namespace aux {

template< typename Predicate, typename LastIterator >
struct iter_fold_if_pred
{
    template< typename State, typename Iterator >
    struct apply
    {
        typedef and_<
              not_< is_same<Iterator,LastIterator> >
            , apply1<Predicate,Iterator>
            > type;
    };
};

} // namespace aux

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(State)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(ForwardOp)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(ForwardPredicate)
    , typename BackwardOp = void_
    , typename BackwardPredicate = void_
    >
struct iter_fold_if
{
private:

    typedef typename begin<Sequence>::type first_;
    typedef typename end<Sequence>::type last_;
    typedef typename lambda<ForwardOp>::type forward_op_;
    typedef typename lambda<ForwardPredicate>::type forward_pred_;
    typedef typename lambda<BackwardOp>::type backward_op_;

    typedef typename apply_if<
          is_void_<BackwardPredicate>
        , if_< is_void_<BackwardOp>, always<false_>, always<true_> >
        , lambda<BackwardPredicate>
        >::type backward_pred_;

// cwpro8 doesn't like 'cut-off' type here (use typedef instead)
#if !BOOST_WORKAROUND(__MWERKS__, BOOST_TESTED_AT(0x3003))
    struct result_ :
#else
    typedef
#endif
        aux::iter_fold_if_impl<
          first_
        , State
        , forward_op_
        , aux::iter_fold_if_pred< forward_pred_,last_ >
        , backward_op_
        , backward_pred_
        >
#if !BOOST_WORKAROUND(__MWERKS__, BOOST_TESTED_AT(0x3003))
    { };
#else
    result_;
#endif

public:

    typedef pair<
          typename result_::state
        , typename result_::iterator
        > type;

    BOOST_MPL_AUX_LAMBDA_SUPPORT(
          6
        , iter_fold_if
        , (Sequence,State,ForwardOp,ForwardPredicate,BackwardOp,BackwardPredicate)
        )
};

BOOST_MPL_AUX_VOID_SPEC(4, iter_fold_if)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_ITER_FOLD_IF_HPP_INCLUDED
