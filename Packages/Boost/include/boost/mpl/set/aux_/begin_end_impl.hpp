
#ifndef BOOST_MPL_SET_AUX_BEGIN_END_IMPL_HPP_INCLUDED
#define BOOST_MPL_SET_AUX_BEGIN_END_IMPL_HPP_INCLUDED

// + file: boost/mpl/aux_/begin_end_impl.hpp
// + last modified: 03/may/03

// Copyright (c) 2002-03
// David Abrahams, Aleksey Gurtovoy
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

#include "boost/mpl/begin_end_fwd.hpp"
#include "boost/mpl/set/aux_/iterator.hpp"

namespace boost { namespace mpl {

template<>
struct begin_traits< aux::set_tag >
{
    template< typename Set > struct algorithm
    {
        typedef set_iter<Set,Set> type;
    };
};

template<>
struct end_traits< aux::set_tag >
{
    template< typename Set > struct algorithm
    {
        typedef set_iter< Set,set0<> > type;
    };
};

}}

#endif // BOOST_MPL_SET_AUX_BEGIN_END_IMPL_HPP_INCLUDED
