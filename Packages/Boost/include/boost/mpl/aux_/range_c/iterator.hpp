//-----------------------------------------------------------------------------
// boost mpl/aux_/range_c/iterator.hpp header file
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

#ifndef BOOST_MPL_AUX_RANGE_C_ITERATOR_HPP_INCLUDED
#define BOOST_MPL_AUX_RANGE_C_ITERATOR_HPP_INCLUDED

#include "boost/mpl/iterator_tag.hpp"
#include "boost/mpl/plus.hpp"
#include "boost/mpl/minus.hpp"
#include "boost/mpl/aux_/iterator_names.hpp"

namespace boost {
namespace mpl {

template< typename N >
struct range_c_iterator
{
    typedef ra_iter_tag_ category;
    typedef N type;

    typedef range_c_iterator<typename N::next> next;
    typedef range_c_iterator<typename N::prior> prior;

    template< typename D >
    struct BOOST_MPL_AUX_ITERATOR_ADVANCE
    {
        typedef range_c_iterator<
              typename plus<N,D>::type
            > type;
    };

    template< typename U >
    struct BOOST_MPL_AUX_ITERATOR_DISTANCE
    {
        typedef typename minus<typename U::type,N>::type type;
    };
};

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_RANGE_C_ITERATOR_HPP_INCLUDED
