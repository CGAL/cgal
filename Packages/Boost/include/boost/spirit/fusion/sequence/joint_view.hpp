/*=============================================================================
    Copyright (c) 2003 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_SEQUENCE_JOINT_VIEW_HPP)
#define FUSION_SEQUENCE_JOINT_VIEW_HPP

#include <boost/spirit/fusion/detail/access.hpp>
#include <boost/spirit/fusion/sequence/begin.hpp>
#include <boost/spirit/fusion/sequence/end.hpp>
#include <boost/spirit/fusion/iterator/joint_view_iterator.hpp>
#include <boost/spirit/fusion/sequence/detail/joint_view_begin_end_traits.hpp>
#include <boost/spirit/fusion/sequence/detail/sequence_base.hpp>
#include <boost/spirit/fusion/sequence/as_fusion_sequence.hpp>

namespace boost { namespace fusion
{
    struct joint_view_tag;

    template <typename View1, typename View2>
    struct joint_view : sequence_base<joint_view<View1, View2> >
    {
        typedef as_fusion_sequence<View1> view1_converter;
        typedef typename view1_converter::type view1;
        typedef as_fusion_sequence<View2> view2_converter;
        typedef typename view2_converter::type view2;

        typedef joint_view_tag tag;
        typedef typename meta::begin<view1>::type first_type;
        typedef typename meta::end<view1>::type last_type;
        typedef typename meta::begin<view2>::type concat_type;
        typedef typename meta::end<view2>::type concat_last_type;

        joint_view(View1& view1, View2& view2)
            : first(fusion::begin(view1))
            , concat(fusion::begin(view2))
            , concat_last(fusion::end(view2))
        {}

        first_type first;
        concat_type concat;
        concat_last_type concat_last;
    };
}}

#endif


