/*=============================================================================
    Copyright (c) 2003 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_ITERATOR_TUPLE_ITERATOR_HPP)
#define FUSION_ITERATOR_TUPLE_ITERATOR_HPP

#include <boost/mpl/int.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/spirit/fusion/iterator/detail/iterator_base.hpp>
#include <boost/spirit/fusion/iterator/detail/tuple_iterator/deref_traits.hpp>
#include <boost/spirit/fusion/iterator/detail/tuple_iterator/value_traits.hpp>
#include <boost/spirit/fusion/iterator/detail/tuple_iterator/next_traits.hpp>
#include <boost/spirit/fusion/iterator/detail/tuple_iterator/prior_traits.hpp>
#include <boost/spirit/fusion/iterator/detail/tuple_iterator/equal_to_traits.hpp>
#include <boost/spirit/fusion/sequence/detail/tuple_begin_end_traits.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/identity.hpp>

namespace boost { namespace fusion
{
    struct tuple_iterator_tag;
    struct void_t;

    template <int N, typename Tuple>
    struct tuple_iterator : iterator_base<tuple_iterator<N, Tuple> >
    {
        typedef mpl::int_<N> index;
        typedef Tuple tuple;
        typedef tuple_iterator_tag tag;
        typedef tuple_iterator<N, Tuple> self_type;

        typedef typename
            mpl::eval_if<
                mpl::less<index, typename Tuple::size>
              , detail::tuple_iterator_next_traits_impl<self_type>
              , mpl::identity<void_t>
            >::type
        next;

        typedef typename
            mpl::eval_if<
                mpl::less<index, typename Tuple::size>
              , detail::tuple_iterator_value_traits_impl<self_type>
              , mpl::identity<void_t>
            >::type
        type;

        tuple_iterator(tuple& t)
            : t(t) {}

        template <int N2, typename Tuple2>
        tuple_iterator(tuple_iterator<N2, Tuple2> const& i)
            : t(static_cast<tuple&>(i.get_tuple())) {}

        tuple&
        get_tuple() const
        {
            return t;
        }

    private:

        tuple& t;
    };
}}

#endif
