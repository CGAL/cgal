//-----------------------------------------------------------------------------
// boost mpl/aux_/single_element_iter.hpp header file
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

#ifndef BOOST_MPL_AUX_SINGLE_ELEMENT_ITER_HPP_INCLUDED
#define BOOST_MPL_AUX_SINGLE_ELEMENT_ITER_HPP_INCLUDED

#include "boost/mpl/iterator_tag.hpp"
#include "boost/mpl/plus.hpp"
#include "boost/mpl/minus.hpp"
#include "boost/mpl/int.hpp"
#include "boost/mpl/aux_/value_wknd.hpp"
#include "boost/mpl/aux_/iterator_names.hpp"
#include "boost/mpl/aux_/lambda_spec.hpp"
#include "boost/mpl/aux_/config/ctps.hpp"
#include "boost/mpl/aux_/config/nttp.hpp"

namespace boost { namespace mpl { 

namespace aux {

template< typename T, int N >
struct single_element_iter;

// random access support
template< typename T, BOOST_MPL_AUX_NTTP_DECL(int, N) >
struct single_iter_base
{
    typedef ra_iter_tag_ category;
    typedef int_<N> position;

    template< typename D >
    struct BOOST_MPL_AUX_ITERATOR_ADVANCE
    {
        typedef plus< int_<N>,D > n_;
        typedef single_element_iter<
              T
            , BOOST_MPL_AUX_VALUE_WKND(n_)::value
            > type;
    };

    template< typename U >
    struct BOOST_MPL_AUX_ITERATOR_DISTANCE
    {
        typedef typename minus<
              typename U::position
            , int_<N>
            >::type type;
    };
};

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

template< typename T >
struct single_element_iter<T,0>
    : single_iter_base<T,0>
{
    typedef single_element_iter<T,1> next;
    typedef T type;
};

template< typename T >
struct single_element_iter<T,1>
    : single_iter_base<T,1>
{
    typedef single_element_iter<T,0> prior;
};

#else

template< BOOST_MPL_AUX_NTTP_DECL(int, N) > struct single_iter_impl
{
    template< typename T > struct result_;
};

template<>
struct single_iter_impl<0>
{
    template< typename T > struct result_
        : single_iter_base<T,0>
    {
        typedef single_element_iter<T,1> next;
        typedef T type;
    };
};

template<>
struct single_iter_impl<1>
{
    template< typename T > struct result_
        : single_iter_base<T,1>
    {
        typedef single_element_iter<T,0> prior;
    };
};

template< typename T, BOOST_MPL_AUX_NTTP_DECL(int, N) >
struct single_element_iter
    : single_iter_impl<N>::template result_<T>
{
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

} // namespace aux

//BOOST_MPL_AUX_PASS_THROUGH_LAMBDA_SPEC(1, aux::single_element_iter)

}} // namespace boost::mpl

#endif // BOOST_MPL_AUX_SINGLE_ELEMENT_ITER_HPP_INCLUDED
