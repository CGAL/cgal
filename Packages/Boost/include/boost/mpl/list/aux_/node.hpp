//-----------------------------------------------------------------------------
// boost mpl/list/aux_/node.hpp header file
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

#ifndef BOOST_MPL_LIST_AUX_NODE_HPP_INCLUDED
#define BOOST_MPL_LIST_AUX_NODE_HPP_INCLUDED

#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/list/aux_/tag.hpp"
#include "boost/mpl/aux_/config/msvc.hpp"
#include "boost/mpl/aux_/config/workaround.hpp"

namespace boost {
namespace mpl {

template<
      typename Size
    , typename T
    , typename Next
    >
struct list_node
{
// agurt, 17/jul/03: to facilitate the deficient 'is_sequence' implementation 
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1300)
    typedef int begin;
#endif
    typedef aux::list_tag tag;
    typedef list_node type;

    typedef Size size;
    typedef T item;
    typedef Next next;
};

struct null_node
{
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1300)
    typedef int begin;
#endif
    typedef aux::list_tag tag;
    typedef null_node type;
    typedef integral_c<long,0> size;
};

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_LIST_AUX_NODE_HPP_INCLUDED
