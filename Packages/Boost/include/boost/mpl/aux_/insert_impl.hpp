//-----------------------------------------------------------------------------
// boost mpl/aux_/insert_impl.hpp header file
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

#ifndef BOOST_MPL_INSERT_IMPL_HPP_INCLUDED
#define BOOST_MPL_INSERT_IMPL_HPP_INCLUDED

#include "boost/mpl/copy_backward.hpp"
#include "boost/mpl/iterator_range.hpp"
#include "boost/mpl/clear.hpp"
#include "boost/mpl/push_front.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/traits_lambda_spec.hpp"
#include "boost/type_traits/is_same.hpp"

namespace boost {
namespace mpl {

// default implementation; conrete sequences might override it by 
// specializing either the |insert_traits| or the primary |insert| template

template< typename Tag >
struct insert_traits
{
    template<
          typename Sequence
        , typename Pos
        , typename T
        >
    struct algorithm
    {
        typedef iterator_range<
              typename begin<Sequence>::type
            , Pos
            > first_half_;

        typedef iterator_range<
              Pos
            , typename end<Sequence>::type
            > second_half_;

        typedef typename copy_backward<
              second_half_
            , typename clear<Sequence>::type
            , push_front<_,_>
            >::type half_sequence_;

        typedef typename copy_backward<
              first_half_
            , typename push_front<half_sequence_,T>::type
            , push_front<_,_>
            >::type type;
    };
};

BOOST_MPL_ALGORITM_TRAITS_LAMBDA_SPEC(3,insert_traits)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_INSERT_IMPL_HPP_INCLUDED
