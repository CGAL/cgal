//-----------------------------------------------------------------------------
// boost mpl/list/aux_/push_front.hpp header file
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

#ifndef BOOST_MPL_LIST_AUX_PUSH_FRONT_HPP_INCLUDED
#define BOOST_MPL_LIST_AUX_PUSH_FRONT_HPP_INCLUDED

#include "boost/mpl/push_front_fwd.hpp"
#include "boost/mpl/aux_/next.hpp"
#include "boost/mpl/list/aux_/node.hpp"
#include "boost/mpl/list/aux_/tag.hpp"

namespace boost {
namespace mpl {

template<>
struct push_front_traits< aux::list_tag >
{
    template< typename List, typename T > struct algorithm
    {
        typedef list_node<
              typename BOOST_MPL_AUX_NEXT(List::size)
            , T
            , typename List::type
            > type;
    };
};

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_LIST_AUX_PUSH_FRONT_HPP_INCLUDED
