//-----------------------------------------------------------------------------
// boost mpl/list/aux_/begin_end.hpp header file
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

#ifndef BOOST_MPL_LIST_AUX_BEGIN_END_HPP_INCLUDED
#define BOOST_MPL_LIST_AUX_BEGIN_END_HPP_INCLUDED

#include "boost/mpl/begin_end_fwd.hpp"
#include "boost/mpl/list/aux_/iterator.hpp"
#include "boost/mpl/list/aux_/tag.hpp"
#include "boost/mpl/list/aux_/node.hpp"

namespace boost {
namespace mpl {

template<>
struct begin_traits< aux::list_tag >
{
    template< typename List > struct algorithm
    {
        typedef list_iterator<typename List::type> type;
    };
};

template<>
struct end_traits< aux::list_tag >
{
    template< typename > struct algorithm
    {
        typedef list_iterator<null_node> type;
    };
};

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_LIST_AUX_BEGIN_END_HPP_INCLUDED
