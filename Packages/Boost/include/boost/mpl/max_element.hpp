//-----------------------------------------------------------------------------
// boost mpl/max_element.hpp header file
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

#ifndef BOOST_MPL_MAX_ELEMENT_HPP_INCLUDED
#define BOOST_MPL_MAX_ELEMENT_HPP_INCLUDED

#include "boost/mpl/less.hpp"
#include "boost/mpl/iter_fold.hpp"
#include "boost/mpl/begin_end.hpp"
#include "boost/mpl/if.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/aux_/common_name_wknd.hpp"
#include "boost/mpl/aux_/void_spec.hpp"

namespace boost {
namespace mpl {

BOOST_MPL_AUX_COMMON_NAME_WKND(max_element)

namespace aux {

template< typename Predicate >
struct select_max
{
    template< typename OldIterator, typename Iterator >
    struct apply
    {
        typedef typename apply2<
              Predicate
            , typename OldIterator::type
            , typename Iterator::type
            >::type condition_;

        typedef typename if_<
              condition_
            , Iterator
            , OldIterator
            >::type type;
    };
};

} // namespace aux 


BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename Predicate = less<_,_>
    >
struct max_element
{
 private:
    typedef typename lambda<Predicate>::type pred_;

 public:
    typedef typename iter_fold<
          Sequence
        , typename begin<Sequence>::type
        , protect< aux::select_max<pred_> >
        >::type type;
};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(1, max_element)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_MAX_ELEMENT_HPP_INCLUDED
