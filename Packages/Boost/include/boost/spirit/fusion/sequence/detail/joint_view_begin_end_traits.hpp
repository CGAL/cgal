/*=============================================================================
    Copyright (c) 2003 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_SEQUENCE_DETAIL_JOINT_VIEW_BEGIN_END_TRAITS_HPP)
#define FUSION_SEQUENCE_DETAIL_JOINT_VIEW_BEGIN_END_TRAITS_HPP

#include <boost/spirit/fusion/detail/config.hpp>
#include <boost/spirit/fusion/iterator/equal_to.hpp>
#include <boost/mpl/if.hpp>

namespace boost { namespace fusion
{
    struct joint_view_tag;

    template <typename View1, typename View2>
    struct joint_view;

    template <typename First, typename Last, typename Concat>
    struct joint_view_iterator;

    namespace meta
    {
        template <typename Tag>
        struct begin_impl;

        template <>
        struct begin_impl<joint_view_tag>
        {
            template <typename Sequence>
            struct apply
            {
                typedef typename Sequence::first_type first_type;
                typedef typename Sequence::last_type last_type;
                typedef typename Sequence::concat_type concat_type;
                typedef meta::equal_to<first_type, last_type> equal_to;

                typedef typename
                    mpl::if_<
                        equal_to
                      , concat_type
                      , joint_view_iterator<first_type, last_type, concat_type>
                    >::type
                type;

                static type
                call(Sequence& s, mpl::true_)
                {
                    return s.concat;
                }

                static type
                call(Sequence& s, mpl::false_)
                {
                    return type(s.first, s.concat);
                }

                static type
                call(Sequence& s)
                {
                    return call(s, equal_to());
                }
            };
        };

        template <typename Tag>
        struct end_impl;

        template <>
        struct end_impl<joint_view_tag>
        {
            template <typename Sequence>
            struct apply
            {
                typedef typename Sequence::concat_last_type type;

                static type
                call(Sequence& s)
                {
                    return s.concat_last;
                }
            };
        };
    }
}}

namespace boost { namespace mpl
{
    template <typename Tag>
    struct begin_impl;

    template <typename Tag>
    struct end_impl;

    template <>
    struct begin_impl<fusion::joint_view_tag>
        : fusion::meta::begin_impl<fusion::joint_view_tag> {};

    template <>
    struct end_impl<fusion::joint_view_tag>
        : fusion::meta::end_impl<fusion::joint_view_tag> {};
}}

#endif


