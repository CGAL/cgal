/*=============================================================================
    Copyright (c) 2003 Joel de Guzman
    Copyright (c) 2004 Peder Holt
    Copyright (c) 2005 Eric Niebler

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_ITERATOR_CONS_ITERATOR_HPP)
#define FUSION_ITERATOR_CONS_ITERATOR_HPP

#include <boost/spirit/fusion/iterator/as_fusion_iterator.hpp>
#include <boost/spirit/fusion/iterator/detail/iterator_base.hpp>
#include <boost/spirit/fusion/sequence/detail/sequence_base.hpp>
#include <boost/spirit/fusion/iterator/detail/cons_iterator/deref_traits.hpp>
#include <boost/spirit/fusion/iterator/detail/cons_iterator/next_traits.hpp>
#include <boost/spirit/fusion/iterator/detail/cons_iterator/value_traits.hpp>

namespace boost { namespace fusion
{
    struct nil;

    struct cons_iterator_tag;

    template <typename Cons = nil>
    struct cons_iterator : iterator_base<cons_iterator<Cons> >
    {
        typedef cons_iterator_tag tag;
        typedef Cons cons_type;
        
        explicit cons_iterator(cons_type& cons_)
            : cons(cons_) {}

        cons_type& cons;
    };

    template <>
    struct cons_iterator<nil> : iterator_base<cons_iterator<nil> >
    {
        typedef cons_iterator_tag tag;
        typedef nil cons_type;
        cons_iterator() {}
        explicit cons_iterator(nil const&) {}
    };

    template <>
    struct cons_iterator<nil const> : iterator_base<cons_iterator<nil const> >
    {
        typedef cons_iterator_tag tag;
        typedef nil const cons_type;
        cons_iterator() {}
        explicit cons_iterator(nil const&) {}
    };
}}

#endif
