//-----------------------------------------------------------------------------
// boost mpl/aux_/joint_iter.hpp header file
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

#ifndef BOOST_MPL_AUX_JOINT_ITER_HPP_INCLUDED
#define BOOST_MPL_AUX_JOINT_ITER_HPP_INCLUDED

#include "boost/mpl/aux_/lambda_spec.hpp"
#include "boost/mpl/aux_/config/ctps.hpp"
#include "boost/type_traits/is_same.hpp"

namespace boost { namespace mpl {

namespace aux {

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

template<
      typename Iterator1
    , typename LastIterator1
    , typename Iterator2
    >
struct joint_iter
{
    typedef Iterator1 base;
    typedef typename base::category category;
    typedef joint_iter<
          typename base::next
        , LastIterator1
        , Iterator2
        > next;

    typedef typename base::type type;
};

template<
      typename LastIterator1
    , typename Iterator2
    >
struct joint_iter<LastIterator1,LastIterator1,Iterator2>
{
    typedef Iterator2 base;
    typedef typename base::category category;
    typedef joint_iter<
          LastIterator1
        , LastIterator1
        , typename base::next
        > next;

    typedef typename base::type type;
};


#else

// forward decl. for 'joint_iter_impl'
template<
      typename Iterator1
    , typename LastIterator1
    , typename Iterator2
    >
struct joint_iter;

template< bool > struct joint_iter_impl
{
    template< 
          typename Iterator1
        , typename LastIterator1
        , typename Iterator2
        >
    struct result_
    {
        typedef Iterator1 base;
        // agurt, 16/nov/02: have to use 'Iterator1' instead of 'base' below
        // to prevent 'base' and 'mpl::base' conflict on MSVC 6.0
        typedef typename Iterator1::category category;
        typedef joint_iter<
              typename Iterator1::next
            , LastIterator1
            , Iterator2
            > next;

        typedef typename Iterator1::type type;
    };
};

template<> struct joint_iter_impl<true>
{
    template< 
          typename Iterator1
        , typename LastIterator1
        , typename Iterator2
        >
    struct result_
    {
        typedef Iterator2 base;
        typedef typename Iterator2::category category;
        typedef joint_iter<
              LastIterator1
            , LastIterator1
            , typename Iterator2::next
            > next;

        typedef typename Iterator2::type type;
    };
};


template<
      typename Iterator1
    , typename LastIterator1
    , typename Iterator2
    >
struct joint_iter
    : joint_iter_impl< is_same<Iterator1,LastIterator1>::value >
        ::template result_<Iterator1,LastIterator1,Iterator2>
{
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

} // namespace aux

BOOST_MPL_AUX_PASS_THROUGH_LAMBDA_SPEC(3, aux::joint_iter)

}} // namespace boost::mpl

#endif // BOOST_MPL_AUX_JOINT_ITER_HPP_INCLUDED
