//-----------------------------------------------------------------------------
// boost mpl/unique.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-02
// Aleksey Gurtovoy and John R. Bandela
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_UNIQUE_HPP_INCLUDED
#define BOOST_MPL_UNIQUE_HPP_INCLUDED

#include "boost/mpl/fold_backward.hpp"
#include "boost/mpl/push_front.hpp"
#include "boost/mpl/clear.hpp"
#include "boost/mpl/front.hpp"
#include "boost/mpl/identity.hpp"
#include "boost/mpl/pair.hpp"
#include "boost/mpl/select1st.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/apply_if.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/void.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_spec.hpp"
#include "boost/mpl/aux_/config/eti.hpp"
#include "boost/type_traits/is_same.hpp"

namespace boost {
namespace mpl {

namespace aux {

template< typename Predicate >
struct unique_op
{
    template< typename Pair, typename T > struct apply
    {
        typedef typename Pair::first seq_;
        typedef typename Pair::second prior_;
        typedef typename apply_if<
              typename apply2<Predicate,prior_,T>::type
            , identity<seq_>
            , push_front<seq_,T>
            >::type new_seq_;

        typedef pair<new_seq_,T> type;
    };
};

} // namespace aux

BOOST_MPL_AUX_PASS_THROUGH_LAMBDA_SPEC(1,aux::unique_op)

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename Predicate = is_same<_,_>
    >
struct unique
{
 private:
    struct none_;
    typedef typename lambda<Predicate>::type pred_;
    typedef typename clear<Sequence>::type result_;
    typedef typename fold_backward<
          Sequence
        , pair<result_,none_>
        , aux::unique_op<pred_>
        >::type fold_result_;

 public:
#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
    // MSVC6.5 forces us to use 'select1st<fold_result_>::type' instead of 
    // simple 'fold_result_::first' here
    typedef typename select1st<fold_result_>::type type;
#else
    typedef typename fold_result_::first type;
#endif
};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(1, unique)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_UNIQUE_HPP_INCLUDED
