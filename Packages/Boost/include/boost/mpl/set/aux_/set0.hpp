
#ifndef BOOST_MPL_SET_AUX_SET0_HPP_INCLUDED
#define BOOST_MPL_SET_AUX_SET0_HPP_INCLUDED

// + file: boost/mpl/aux_/set0.hpp
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

#include "boost/mpl/void_fwd.hpp"
#include "boost/mpl/set/aux_/tag.hpp"
#include "boost/mpl/aux_/config/static_constant.hpp"

namespace boost {
namespace mpl {

template< int p_ = 0 > struct set0
{
    typedef aux::set_tag    tag;
    typedef void_           last_masked;
    typedef void_           item;

    BOOST_STATIC_CONSTANT(long, order = 0);
};

}}

#endif // BOOST_MPL_SET_AUX_SET0_HPP_INCLUDED
