//-----------------------------------------------------------------------------
// boost mpl/count_if.hpp header file
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

#ifndef BOOST_MPL_COUNT_IF_HPP_INCLUDED
#define BOOST_MPL_COUNT_IF_HPP_INCLUDED

#include "boost/mpl/fold.hpp"
#include "boost/mpl/next.hpp"
#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/identity.hpp"
#include "boost/mpl/apply_if.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/protect.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"

namespace boost {
namespace mpl {

namespace aux {

template< typename Predicate >
struct next_if
{
    template<
          typename N
        , typename T
        >
    struct apply
    {
        typedef typename apply_if<
              typename apply1<Predicate,T>::type
            , next<N>
            , identity<N>
            >::type type;
    };
};

} // namespace aux

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Predicate)
    >
struct count_if
{
 private:
    typedef typename lambda<Predicate>::type pred_;

 public:
    typedef typename fold<
          Sequence
        , integral_c<unsigned long, 0>
        , protect< aux::next_if<pred_> >
        >::type type;

    BOOST_MPL_AUX_LAMBDA_SUPPORT(2,count_if,(Sequence,Predicate))
};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(2, count_if)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_COUNT_IF_HPP_INCLUDED
