/*=============================================================================
    Copyright (c) 2003 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_SEQUENCE_AS_FUSION_SEQUENCE_HPP)
#define FUSION_SEQUENCE_AS_FUSION_SEQUENCE_HPP

#include <boost/spirit/fusion/sequence/is_sequence.hpp>
#include <boost/spirit/fusion/sequence/type_sequence.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/static_assert.hpp>

namespace boost { namespace fusion
{
    //  Test T. If it is a fusion sequence, return a reference to it.
    //  else, assume it is an mpl sequence. Fail if it is not.

    template <typename T>
    struct as_fusion_sequence
    {
        typedef typename
            mpl::if_<
                fusion::is_sequence<T>
              , T
              , type_sequence<T>
            >::type
        type;

        static T&
        convert(T& x, mpl::true_)
        {
            return x;
        }

        static type_sequence<T>
        convert(T& x, mpl::false_)
        {
            BOOST_STATIC_ASSERT(mpl::is_sequence<T>::value);
            return type_sequence<T>();
        }

        static typename
            mpl::if_<
                fusion::is_sequence<T>
              , T&
              , type_sequence<T>
            >::type
        convert(T& x)
        {
            return convert(x, fusion::is_sequence<T>());
        }

        static T const&
        convert_const(T const& x, mpl::true_)
        {
            return x;
        }

        static type_sequence<T>
        convert_const(T const& x, mpl::false_)
        {
            BOOST_STATIC_ASSERT(mpl::is_sequence<T>::value);
            return type_sequence<T>();
        }

        static typename
            mpl::if_<
                fusion::is_sequence<T>
              , T const&
              , type_sequence<T>
            >::type
        convert_const(T const& x)
        {
            return convert_const(x, fusion::is_sequence<T>());
        }
    };
}}

#endif
