/*=============================================================================
    Copyright (c) 2003 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_ADAPT_DEREF_TRAITS_HPP)
#define FUSION_ADAPT_DEREF_TRAITS_HPP

#include <boost/spirit/fusion/iterator/deref.hpp>

namespace boost { namespace fusion { namespace detail
{
    struct adapt_deref_traits
    {
        template <typename Iterator>
        struct apply
        {
            typedef typename
                meta::deref<typename Iterator::first_type>::type
            type;

            static type
            call(Iterator const& i)
            {
                return *i.first;
            }
        };
    };
}}}

#endif


