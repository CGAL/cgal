//-----------------------------------------------------------------------------
// boost mpl/list/list0_c.hpp header file
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

#ifndef BOOST_MPL_LIST_LIST0_C_HPP_INCLUDED
#define BOOST_MPL_LIST_LIST0_C_HPP_INCLUDED

#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/list/aux_/push_front.hpp"
#include "boost/mpl/list/aux_/pop_front.hpp"
#include "boost/mpl/list/aux_/front.hpp"
#include "boost/mpl/list/aux_/clear.hpp"
#include "boost/mpl/list/aux_/size.hpp"
#include "boost/mpl/list/aux_/empty.hpp"
#include "boost/mpl/list/aux_/begin_end.hpp"
#include "boost/mpl/list/aux_/node.hpp"

namespace boost {
namespace mpl {

template< typename T > struct list0_c
    : null_node
{
    typedef null_node type;
};

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_LIST_LIST0_C_HPP_INCLUDED
