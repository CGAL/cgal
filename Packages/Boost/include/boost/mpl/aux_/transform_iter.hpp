//-----------------------------------------------------------------------------
// boost mpl/aux_/transform_iter.hpp header file
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

#ifndef BOOST_MPL_AUX_TRANSFORM_ITER_HPP_INCLUDED
#define BOOST_MPL_AUX_TRANSFORM_ITER_HPP_INCLUDED

#include "boost/mpl/apply.hpp"
#include "boost/mpl/aux_/lambda_spec.hpp"
#include "boost/mpl/aux_/config/ctps.hpp"
#include "boost/type_traits/is_same.hpp"

namespace boost {
namespace mpl {
namespace aux {

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

template<
      typename Iterator
    , typename LastIterator
    , typename F
    >
struct transform_iter
{
    typedef Iterator base;
    typedef typename base::category category;
    typedef transform_iter< typename base::next,LastIterator,F > next;
    
    typedef typename apply1<
          F
        , typename base::type
        >::type type;
};

template<
      typename LastIterator
    , typename F
    >
struct transform_iter< LastIterator,LastIterator,F >
{
    typedef LastIterator base;
    typedef typename base::category category;
};

#else

template<
      typename Iterator
    , typename LastIterator
    , typename F
    >
struct transform_iter;

template< bool >
struct transform_iter_impl 
{
    template<
          typename Iterator
        , typename LastIterator
        , typename F
        >
    struct result_
    {
        typedef Iterator base;
        // agurt, 14/oct/02: have to use |Iterator| instead of |base| below
        // to prevent |base| and |mpl::base| conflict on MSVC 6.0
        typedef typename Iterator::category category;
        typedef transform_iter< typename Iterator::next,LastIterator,F > next;
        
        typedef typename apply1<
              F
            , typename Iterator::type
            >::type type;
    };
};

template<>
struct transform_iter_impl<true>
{
    template<
          typename Iterator
        , typename LastIterator
        , typename F
        >
    struct result_
    {
        typedef Iterator base;
        typedef typename Iterator::category category;
    };
};

template<
      typename Iterator
    , typename LastIterator
    , typename F
    >
struct transform_iter
    : transform_iter_impl<
          ::boost::is_same<Iterator,LastIterator>::value
        >::template result_< Iterator,LastIterator,F >
{
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

} // namespace aux

BOOST_MPL_AUX_PASS_THROUGH_LAMBDA_SPEC(3, aux::transform_iter)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_AUX_TRANSFORM_ITER_HPP_INCLUDED
