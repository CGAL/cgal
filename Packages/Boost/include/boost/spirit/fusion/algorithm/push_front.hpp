/*=============================================================================
    Copyright (c) 2003 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_ALGORITHM_PUSH_FRONT_HPP)
#define FUSION_ALGORITHM_PUSH_FRONT_HPP

#include <boost/spirit/fusion/sequence/single_view.hpp>
#include <boost/spirit/fusion/sequence/joint_view.hpp>
#include <boost/spirit/fusion/sequence/detail/sequence_base.hpp>

namespace boost { namespace fusion
{
    namespace meta
    {
        template <typename Sequence, typename T>
        struct push_front
        {
            typedef joint_view<single_view<T>, Sequence> type;
        };
    }

    namespace function
    {
        struct push_front
        {
            template <typename Sequence, typename T>
            struct apply : meta::push_front<Sequence, T> {};

            template <typename Sequence, typename T>
            inline typename apply<Sequence const, T>::type
            operator()(Sequence const& seq, T const& x) const
            {
                typedef joint_view<single_view<T>, Sequence const> result;
                single_view<T> val(x);
                return result(val, seq);
            }

            template <typename Sequence, typename T>
            inline typename apply<Sequence, T>::type
            operator()(Sequence& seq, T const& x) const
            {
                typedef joint_view<single_view<T>, Sequence> result;
                single_view<T> val(x);
                return result(val, seq);
            }
        };
    }

    function::push_front const push_front = function::push_front();
}}

#endif

