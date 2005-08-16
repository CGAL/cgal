
//  (C) Copyright Daniel James 2005.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Based on Peter Dimov's proposal
//  http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2005/n1756.pdf
//  issue 6.18. 

#if !defined(BOOST_FUNCTIONAL_HASH_SET_HPP)
#define BOOST_FUNCTIONAL_HASH_SET_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <boost/config.hpp>
#include <set>
#include <boost/functional/hash/hash.hpp>

namespace boost
{
    template <class K, class C, class A>
    std::size_t hash_value(std::set<K, C, A> const& v)
    {
        return hash_range(v.begin(), v.end());
    }

    template <class K, class C, class A>
    std::size_t hash_value(std::multiset<K, C, A> const& v)
    {
        return hash_range(v.begin(), v.end());
    }

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
    namespace hash_detail
    {
        template <class K, class C, class A>
        struct call_hash<std::set<K, C, A> >
        {
            static std::size_t call(std::set<K, C, A> const& val)
            {
                return boost::hash_value(val);
            }
        };

        template <class K, class C, class A>
        struct call_hash<std::multiset<K, C, A> >
        {
            static std::size_t call(std::multiset<K, C, A> const& val)
            {
                return boost::hash_value(val);
            }
        };
    }
#endif
}

#endif

