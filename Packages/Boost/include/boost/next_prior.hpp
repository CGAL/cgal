//  Boost next_prior.hpp header file  ---------------------------------------//

//  (C) Copyright Boost.org 1999-2003. Permission to copy, use, modify, sell
//  and distribute this software is granted provided this copyright
//  notice appears in all copies. This software is provided "as is" without
//  express or implied warranty, and with no claim as to its suitability for
//  any purpose.

//  See http://www.boost.org/libs/utility for documentation.

//  Revision History
//  13 Dec 2003  Added next(x, n) and prior(x, n) (Daniel Walker)

#ifndef BOOST_NEXT_PRIOR_HPP_INCLUDED
#define BOOST_NEXT_PRIOR_HPP_INCLUDED

#include <iterator>

namespace boost {

//  Helper functions for classes like bidirectional iterators not supporting
//  operator+ and operator-
//
//  Usage:
//    const std::list<T>::iterator p = get_some_iterator();
//    const std::list<T>::iterator prev = boost::prior(p);
//    const std::list<T>::iterator next = boost::next(prev, 2);

//  Contributed by Dave Abrahams

template <class T>
inline T next(T x) { return ++x; }

template <class T, class Distance>
inline T next(T x, Distance n)
{
    std::advance(x, n);
    return x;
}

template <class T>
inline T prior(T x) { return --x; }

template <class T, class Distance>
inline T prior(T x, Distance n)
{
    std::advance(x, -n);
    return x;
}

} // namespace boost

#endif  // BOOST_NEXT_PRIOR_HPP_INCLUDED
