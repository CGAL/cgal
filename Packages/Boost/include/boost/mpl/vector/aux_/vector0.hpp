//-----------------------------------------------------------------------------
// boost mpl/vector/aux_/vector0.hpp header file
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

#ifndef BOOST_MPL_VECTOR_AUX_VECTOR0_HPP_INCLUDED
#define BOOST_MPL_VECTOR_AUX_VECTOR0_HPP_INCLUDED

#include "boost/mpl/integral_c.hpp"
#include "boost/mpl/void.hpp"
#include "boost/mpl/aux_/type_wrapper.hpp"

#include "boost/mpl/vector/aux_/iterator.hpp"
#include "boost/mpl/vector/aux_/tag.hpp"
#include "boost/mpl/aux_/config/vector.hpp"

namespace boost {
namespace mpl {

#if defined(BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL)

template< typename Dummy = void_ > struct vector0;
template<> struct vector0<void_>
{
    static aux::type_wrapper<void_> item(...);

    typedef aux::vector_tag tag;
    typedef integral_c<long,0> size;
    typedef vector0 type;
};

#else

template< typename Dummy = void_ > struct vector0;
template<> struct vector0<void_>
{
    typedef aux::vector_tag<0> tag;
    typedef vector0 type;
    typedef void_ item0;
    
    typedef vector_iterator< vector0,integral_c<long,0> > begin;
    typedef vector_iterator< vector0,integral_c<long,0> > end;
};

#endif // BOOST_MPL_TYPEOF_BASED_VECTOR_IMPL

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_VECTOR_AUX_VECTOR0_HPP_INCLUDED
