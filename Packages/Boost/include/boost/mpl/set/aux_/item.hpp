
#ifndef BOOST_MPL_SET_AUX_ITEM_HPP_INCLUDED
#define BOOST_MPL_SET_AUX_ITEM_HPP_INCLUDED

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

#include "boost/mpl/long.hpp"
#include "boost/mpl/void_fwd.hpp"
#include "boost/mpl/set/aux_/set0.hpp"
#include "boost/mpl/aux_/yes_no.hpp"
#include "boost/mpl/aux_/type_wrapper.hpp"
#include "boost/mpl/aux_/config/static_constant.hpp"
#include "boost/mpl/aux_/config/workaround.hpp"

namespace boost {
namespace mpl {

aux::no_tag operator/(set0<> const&, void*);
aux::no_tag operator%(set0<> const&, void*);

// agurt, 03/may/03: forward declarations of operators, to supressing a GCC warning, 
// see below; breaks 2.95.x!
#if BOOST_WORKAROUND(__GNUC__, > 2) && BOOST_WORKAROUND(__GNUC__, BOOST_TESTED_AT(3))
template< typename T, typename Base > struct s_item;
template< typename T, typename Base > struct s_mask;

template< typename T, typename Base >
typename s_item<T,Base>::order_tag
operator/(s_item<T,Base> const&, aux::type_wrapper<T>*);

template< typename T, typename Base >
aux::yes_tag operator%(s_item<T,Base> const&, aux::type_wrapper<T>*);

template< typename T, typename Base >
aux::no_tag operator%(s_mask<T,Base> const&, aux::type_wrapper<T>*);
#endif

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x561))
namespace aux {
template< long n_ > struct order_tag_
{
    typedef char (&type)[n_];    
};
}
#endif

template< typename T, typename Base >
struct s_item
    : Base
{
    typedef void_   last_masked;
    typedef T       item;
    typedef Base    base;
    
    BOOST_STATIC_CONSTANT(long, order = Base::order + 1);

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x561))
    typedef typename aux::order_tag_<(Base::order + 1)>::type order_tag;
#else
    typedef char (&order_tag)[order];
#endif

#if BOOST_WORKAROUND(__GNUC__, > 2) && BOOST_WORKAROUND(__GNUC__, BOOST_TESTED_AT(3))
    // to make GCC happy
    friend order_tag    operator/<>(s_item const&, aux::type_wrapper<T>*);
    friend aux::yes_tag operator%<>(s_item const&, aux::type_wrapper<T>*);
#else
    friend order_tag    operator/(s_item const&, aux::type_wrapper<T>*);
    friend aux::yes_tag operator%(s_item const&, aux::type_wrapper<T>*);
#endif
};


template< typename T, typename Base >
struct s_mask
    : Base
{
    typedef T last_masked;
#if BOOST_WORKAROUND(__GNUC__, > 2) && BOOST_WORKAROUND(__GNUC__, BOOST_TESTED_AT(3))
    // to make GCC happy
    friend aux::no_tag operator%<>(s_mask const&, aux::type_wrapper<T>*);
#else
    friend aux::no_tag operator%(s_mask const&, aux::type_wrapper<T>*);
#endif
};

}}

#endif // BOOST_MPL_SET_AUX_ITEM_HPP_INCLUDED
