//-----------------------------------------------------------------------------
// boost mpl/distance.hpp header file
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

#ifndef BOOST_MPL_DISTANCE_HPP_INCLUDED
#define BOOST_MPL_DISTANCE_HPP_INCLUDED

#include "boost/mpl/aux_/iter_distance.hpp"
#include "boost/mpl/aux_/iterator_category.hpp"
#include "boost/mpl/iterator_tag.hpp"
#include "boost/mpl/iter_fold.hpp"
#include "boost/mpl/iterator_range.hpp"
#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/next.hpp"
#include "boost/mpl/aux_/common_name_wknd.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/config.hpp"

namespace boost {
namespace mpl {

BOOST_MPL_AUX_COMMON_NAME_WKND(distance)

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace aux {

// forward/bidirectional iterators
template< typename Category, typename First, typename Last >
struct distance_impl
    : iter_fold<
          iterator_range<First,Last>
        , integral_c<long, 0>
        , next<>
        >
{
};

template< typename First, typename Last >
struct distance_impl<ra_iter_tag_,First,Last>
    : aux::iter_distance<First,Last>
{
};

} // namespace aux

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(First)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Last)
    >
struct distance
// Aleksey claims borland doesn't like inheritance here, but it passes
// all the same tests; the workaround can easily be enabled again if
// verified.  -- dwa 2003/5/8
# if 1 || !defined(__BORLANDC__)  
    : aux::distance_impl<
          typename BOOST_MPL_AUX_ITERATOR_CATEGORY(First)
        , First
        , Last
      >
# endif 
{
# if 0 && defined(__BORLANDC__)
    typedef typename aux::distance_impl<
        typename BOOST_MPL_AUX_ITERATOR_CATEGORY(First)
      , First
      , Last
    >::type type;
    BOOST_STATIC_CONSTANT(typename type::value_type, value = type::value);
# endif 
};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

#else

namespace aux {

// forward/bidirectional iterators
template< typename Category >
struct distance_impl
{
    template< typename First, typename Last > struct result_
        : iter_fold<
              iterator_range<First,Last>
            , integral_c<long, 0>
            , next<>
            >
    {
    };
};

template<>
struct distance_impl<ra_iter_tag_>
{
    template< typename First, typename Last > struct result_
        : aux::iter_distance<First,Last>
    {
    };
};

} // namespace aux

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(First)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Last)
    >
struct distance
#if !defined(BOOST_MSVC) || BOOST_MSVC != 1300
    : aux::distance_impl< typename BOOST_MPL_AUX_ITERATOR_CATEGORY(First) >
#else
    : aux::distance_impl< fwd_iter_tag_ >
#endif
        ::template result_<First,Last>
{
};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(2, distance)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_DISTANCE_HPP_INCLUDED
