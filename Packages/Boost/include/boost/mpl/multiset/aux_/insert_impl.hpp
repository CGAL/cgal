
#ifndef BOOST_MPL_MULTISET_AUX_INSERT_IMPL_HPP_INCLUDED
#define BOOST_MPL_MULTISET_AUX_INSERT_IMPL_HPP_INCLUDED

// + file: boost/mpl/multiset/aux_/insert_impl.hpp
// + last modified: 05/nov/03

// Copyright Aleksey Gurtovoy 2003
//
// Use, modification and distribution are subject to the Boost Software 
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy 
// at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

#include "boost/mpl/multiset/aux_/item.hpp"
#include "boost/mpl/multiset/aux_/tag.hpp"
#include "boost/mpl/insert_fwd.hpp"

namespace boost { namespace mpl {
//#error here!
template<>
struct insert_traits< aux::multiset_tag >
{
    template< typename Set, typename Key, typename unused_ > struct algorithm
    {
        typedef ms_item<Key,Set> type;
    };
};

}} // namespace boost::mpl

#endif // BOOST_MPL_MULTISET_AUX_INSERT_IMPL_HPP_INCLUDED
