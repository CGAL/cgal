//-----------------------------------------------------------------------------
// boost mpl/transform.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-03
// Dave Abrahams, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_TRANSFORM_HPP_INCLUDED
#define BOOST_MPL_TRANSFORM_HPP_INCLUDED

#include "boost/mpl/fold_backward.hpp"
#include "boost/mpl/push_front.hpp"
#include "boost/mpl/clear.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/pair.hpp"
#include "boost/mpl/protect.hpp"
#include "boost/mpl/iterator_tag.hpp"
#include "boost/mpl/apply_if.hpp"
#include "boost/mpl/aux_/common_name_wknd.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/type_traits/is_same.hpp"

namespace boost {
namespace mpl {

BOOST_MPL_AUX_COMMON_NAME_WKND(transform)

namespace aux {

template< typename Op >
struct transform_op
{
    template< typename Sequence, typename T > struct apply
    {
        typedef typename push_front<
              Sequence
            , typename apply1<Op,T>::type
            >::type type;
    };
};

template< typename Op >
struct transform2_op
{
    template< typename Sequence, typename T > struct apply
    {
        typedef typename push_front<
              Sequence
            , typename apply2<
                    Op
                  , typename T::first
                  , typename T::second
                  >::type
            >::type type;
    };
};

template< typename I1, typename I2 >
struct pair_iterator
{
    typedef input_iter_tag_ category;
    typedef pair<typename I1::type, typename I2::type> type;
    typedef pair_iterator<typename I1::next, typename I2::next> next;
};

template<
      typename Seq1, typename Seq2
    >
struct pair_view
{
 public:
    typedef nested_begin_end_tag tag;
    typedef pair_iterator<
        typename mpl::begin<Seq1>::type
      , typename mpl::begin<Seq2>::type
    > begin;
    
    typedef pair_iterator<
        typename mpl::end<Seq1>::type
      , typename mpl::end<Seq2>::type
    > end;
};

} // namespace aux

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Operation)
    >
struct transform1
{
 private:
    typedef typename lambda<Operation>::type op_;
    typedef typename clear<Sequence>::type result_;

 public:
    typedef typename fold_backward<
          Sequence
        , result_
        , protect< aux::transform_op<op_> >
        >::type type;
};

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Seq1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Seq2)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Operation)
    >
struct transform2
{
 private:
    typedef typename lambda<Operation>::type op_;
    typedef typename clear<Seq1>::type result_;

 public:
    typedef typename fold_backward<
          aux::pair_view<Seq1,Seq2>
        , result_
        , protect< aux::transform2_op<op_> >
        >::type type;
};

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Seq1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Seq2OrOperation)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Operation)
    >
struct transform
  : apply_if<
        is_same<Operation,void_>
      , transform1<Seq1,Seq2OrOperation>
      , transform2<Seq1,Seq2OrOperation,Operation>
    >
{
};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(2, transform1)
BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(3, transform2)
BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(3, transform)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_TRANSFORM_HPP_INCLUDED
