
#ifndef BOOST_MPL_SET_AUX_ORDER_IMPL_HPP_INCLUDED
#define BOOST_MPL_SET_AUX_ORDER_IMPL_HPP_INCLUDED

// + file: boost/mpl/aux_/item.hpp
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

#include "boost/mpl/set/aux_/tag.hpp"
#include "boost/mpl/order_fwd.hpp"
#include "boost/mpl/long.hpp"
#include "boost/mpl/void.hpp"
#include "boost/mpl/aux_/static_cast.hpp"
#include "boost/mpl/aux_/type_wrapper.hpp"
#include "boost/mpl/aux_/ptr_to_ref.hpp"
#include "boost/mpl/aux_/config/static_constant.hpp"
#include "boost/mpl/aux_/config/workaround.hpp"

namespace boost {
namespace mpl {

namespace aux {
template< long n_ > struct long_or_void : long_<n_> {};
template<>          struct long_or_void<0> : void_ {};
}

template<>
struct order_impl< aux::set_tag >
{

    template< typename Set, typename T > struct apply
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1300) || BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x561))
        {
        BOOST_STATIC_CONSTANT(long, value = 
              sizeof( 
                  *BOOST_MPL_AUX_STATIC_CAST(Set*, 0)
                    / BOOST_MPL_AUX_STATIC_CAST(aux::type_wrapper<T>*, 0)
                ) - 1
            );

#   if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x561))
        typedef typename aux::long_or_void<(apply::value)>::type type;
#   else
        typedef typename aux::long_or_void<value>::type type;
#   endif

#else
        : aux::long_or_void< 
              sizeof( 
                  aux::ptr_to_ref(BOOST_MPL_AUX_STATIC_CAST(Set*, 0))
                    / BOOST_MPL_AUX_STATIC_CAST(aux::type_wrapper<T>*, 0)
                ) - 1
            >
    {
#endif
    };

};

}}

#endif // BOOST_MPL_SET_AUX_ORDER_IMPL_HPP_INCLUDED
