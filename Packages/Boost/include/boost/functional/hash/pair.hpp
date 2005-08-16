
//  (C) Copyright Daniel James 2005.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Based on Peter Dimov's proposal
//  http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2005/n1756.pdf
//  issue 6.18. 

#if !defined(BOOST_FUNCTIONAL_HASH_PAIR_HPP)
#define BOOST_FUNCTIONAL_HASH_PAIR_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <boost/config.hpp>
#include <utility>
#include <boost/functional/hash/hash.hpp>

namespace boost
{
    template <class A, class B>
    std::size_t hash_value(std::pair<A, B> const& v)
    {
        std::size_t seed = 0;
        hash_combine(seed, v.first);
        hash_combine(seed, v.second);
        return seed;
    }

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
    namespace hash_detail
    {
        template <class A, class B>
        struct call_hash<std::pair<A, B> >
        {
            static std::size_t call(std::pair<A, B> const& val)
            {
                return boost::hash_value(val);
            }
        };
    }
#endif
}

#endif
