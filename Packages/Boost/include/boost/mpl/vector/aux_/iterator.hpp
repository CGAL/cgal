
#ifndef BOOST_MPL_AUX_VECTOR_ITERATOR_HPP_INCLUDED
#define BOOST_MPL_AUX_VECTOR_ITERATOR_HPP_INCLUDED

// + file: boost/mpl/aux_/vector/iterator.hpp
// + last modified: 30/may/03

// Copyright (c) 2000-03
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.
//
// See http://www.boost.org/libs/mpl for documentation.

#include "boost/mpl/iterator_tag.hpp"
#include "boost/mpl/plus.hpp"
#include "boost/mpl/minus.hpp"
#include "boost/mpl/vector/aux_/item.hpp"
#include "boost/mpl/aux_/iterator_names.hpp"
#include "boost/mpl/aux_/value_wknd.hpp"
#include "boost/mpl/aux_/config/workaround.hpp"

namespace boost {
namespace mpl {

template<
      typename Vector
    , typename Pos
    >
struct vector_iterator
{
    typedef ra_iter_tag_ category;
    typedef typename vector_item<
          Vector
        , BOOST_MPL_AUX_VALUE_WKND(Pos)::value
        >::type type;

    typedef Pos pos;
    typedef vector_iterator<Vector,typename Pos::next> next;
    typedef vector_iterator<Vector,typename Pos::prior> prior;

    template< typename Distance >
    struct BOOST_MPL_AUX_ITERATOR_ADVANCE
    {
        typedef vector_iterator<
              Vector
            , typename plus<Pos,Distance>::type
            > type;
    };

    template< typename Other >
    struct BOOST_MPL_AUX_ITERATOR_DISTANCE
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
        : minus<typename Other::pos,Pos>
    {
#else
    {
        typedef typename minus<typename Other::pos,Pos>::type type;
        BOOST_STATIC_CONSTANT(typename Pos::value_type
            , value = (minus<typename Other::pos,Pos>::value));
#endif
    };
};

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_VECTOR_ITERATOR_HPP_INCLUDED
