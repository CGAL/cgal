//-----------------------------------------------------------------------------
// boost mpl/is_sequence.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_IS_SEQUENCE_HPP_INCLUDED
#define BOOST_MPL_IS_SEQUENCE_HPP_INCLUDED

#include "boost/mpl/not.hpp"
#include "boost/mpl/and.hpp"
#include "boost/mpl/begin_end.hpp"
#include "boost/mpl/if.hpp"
#include "boost/mpl/bool.hpp"
#include "boost/mpl/sequence_tag_fwd.hpp"
#include "boost/mpl/identity.hpp"
#include "boost/mpl/void.hpp"
#include "boost/mpl/aux_/has_tag.hpp"
#include "boost/mpl/aux_/has_begin.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"
#include "boost/mpl/aux_/config/eti.hpp"
#include "boost/mpl/aux_/config/workaround.hpp"

#include "boost/type_traits/is_same.hpp"
#include "boost/type_traits/is_class.hpp"

namespace boost { namespace mpl {

#if BOOST_WORKAROUND(BOOST_MSVC, <= 1300)

namespace aux {

// agurt, 11/jun/03: 
// MSVC 6.5/7.0 fails if 'has_begin' is instantiated on a class type that has a
// 'begin' member that doesn't name a type; e.g. 'has_begin< std::vector<int> >'
// would fail; requiring 'T' to have _both_ 'tag' and 'begin' members workarounds
// the issue for most real-world cases
template< typename T > struct is_sequence_impl
    : and_<
          identity< aux::has_tag<T> >
        , identity< aux::has_begin<T> >
        >
{
};

} // namespace aux
        
template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    >
struct is_sequence
    : if_<
          boost::is_class<T> 
        , aux::is_sequence_impl<T>
        , bool_<false>
        >::type
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1,is_sequence,(T))
};

#elif defined(BOOST_MPL_NO_AUX_HAS_XXX)

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    >
struct is_sequence
    : bool_<false>
{
};

#else

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    >
struct is_sequence
    : not_< is_same< typename begin<T>::type, void_ > >
{
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1,is_sequence,(T))
};

#endif // BOOST_MSVC

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
template<> struct is_sequence<int>
    : bool_<false>
{
};
#endif

BOOST_MPL_AUX_VOID_SPEC(1, is_sequence)

}} // namespace boost::mpl

#endif // BOOST_MPL_IS_SEQUENCE_HPP_INCLUDED
