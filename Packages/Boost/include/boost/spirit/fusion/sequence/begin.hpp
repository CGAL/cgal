/*=============================================================================
    Copyright (c) 2003 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_SEQUENCE_BEGIN_HPP)
#define FUSION_SEQUENCE_BEGIN_HPP

#include <boost/spirit/fusion/detail/config.hpp>
#include <boost/spirit/fusion/sequence/as_fusion_sequence.hpp>

namespace boost { namespace fusion
{
    namespace meta
    {
        template <typename Tag>
        struct begin_impl
        {
            template <typename Sequence>
            struct apply;
        };

        template <typename Sequence>
        struct begin
        {
            typedef as_fusion_sequence<Sequence> seq_converter;
            typedef typename seq_converter::type seq;

            typedef typename
                begin_impl<typename seq::tag>::
                    template apply<seq>::type
            type;
        };
    }

    template <typename Sequence>
    inline typename meta::begin<Sequence const>::type
    begin(Sequence const& seq)
    {
        typedef meta::begin<Sequence const> begin_meta;
        return meta::begin_impl<typename begin_meta::seq::tag>::
            template apply<typename begin_meta::seq const>::call(
                begin_meta::seq_converter::convert_const(seq));
    }

    template <typename Sequence>
    inline typename meta::begin<Sequence>::type
    begin(Sequence& seq)
    {
        typedef meta::begin<Sequence> begin_meta;
        return meta::begin_impl<typename begin_meta::seq::tag>::
            template apply<typename begin_meta::seq>::call(
                begin_meta::seq_converter::convert(seq));
    }
}}

#endif
