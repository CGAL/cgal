//-----------------------------------------------------------------------------
// boost mpl/aux_/sort_impl.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002-2003
// Eric Friedman
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_SORT_IMPL_HPP_INCLUDED
#define BOOST_MPL_AUX_SORT_IMPL_HPP_INCLUDED

#include "boost/mpl/aux_/select1st_wknd.hpp"
#include "boost/mpl/aux_/select2nd_wknd.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/apply_if.hpp"
#include "boost/mpl/copy_backward.hpp"
#include "boost/mpl/empty.hpp"
#include "boost/mpl/front.hpp"
#include "boost/mpl/identity.hpp"
#include "boost/mpl/partition.hpp"
#include "boost/mpl/pop_front.hpp"
#include "boost/mpl/push_front.hpp"
#include "boost/mpl/aux_/traits_lambda_spec.hpp"

namespace boost {
namespace mpl {

namespace aux {

template < typename Sequence, typename Predicate > struct quick_sort;

template <typename Predicate, typename Pivot>
struct quick_sort_pred
{
    template <typename T>
    struct apply
    {
        typedef typename apply2< Predicate, T, Pivot >::type
            type;
    };
};

template <typename Sequence, typename Predicate>
struct quick_sort_impl
{
private:

    typedef typename front<Sequence>::type pivot_;
    typedef typename pop_front<Sequence>::type seq_;

    typedef typename partition<
          seq_
        , protect< quick_sort_pred<Predicate,pivot_> >
        >::type partitioned;

    typedef typename quick_sort<
          typename BOOST_MPL_AUX_SELECT1ST_WKND(partitioned), Predicate
        >::type first_part;
    typedef typename quick_sort<
          typename BOOST_MPL_AUX_SELECT2ND_WKND(partitioned), Predicate
        >::type second_part;

public:

    typedef typename copy_backward<
          first_part
        , typename push_front< second_part,pivot_ >::type
        , push_front<_,_>
        >::type type;

};

template <typename Sequence, typename Predicate>
struct quick_sort
    : apply_if<
          empty<Sequence>
        , identity< Sequence >
        , quick_sort_impl< Sequence,Predicate >
        >
{
};

} // namespace aux

template< typename Tag >
struct sort_traits
{
    template< typename Sequence, typename Predicate >
    struct algorithm
    {
        typedef typename aux::quick_sort<
              Sequence, Predicate
            >::type type;
    };
};

BOOST_MPL_ALGORITM_TRAITS_LAMBDA_SPEC(2,sort_traits)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_SORT_IMPL_HPP_INCLUDED
