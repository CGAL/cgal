#ifndef BOOST_DETAIL_ATOMIC_COUNT_GCC_HPP_INCLUDED
#define BOOST_DETAIL_ATOMIC_COUNT_GCC_HPP_INCLUDED

//
//  boost/detail/atomic_count_gcc.hpp
//
//  atomic_count for GNU libstdc++ v3
//
//  http://gcc.gnu.org/onlinedocs/porting/Thread-safety.html
//
//  Copyright (c) 2001, 2002 Peter Dimov and Multi Media Ltd.
//  Copyright (c) 2002 Lars Gullik Bjønnes <larsbj@lyx.org>
//
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//

#include <bits/atomicity.h>

namespace boost
{

namespace detail
{

class atomic_count
{
public:

    explicit atomic_count(long v) : value_(v) {}

    void operator++()
    {
        __atomic_add(&value_, 1);
    }

    long operator--()
    {
        return !__exchange_and_add(&value_, -1);
    }

    operator long() const
    {
        return __exchange_and_add(&value_, 0);
    }

private:

    atomic_count(atomic_count const &);
    atomic_count & operator=(atomic_count const &);

    _Atomic_word value_;
};

} // namespace detail

} // namespace boost

#endif // #ifndef BOOST_DETAIL_ATOMIC_COUNT_GCC_HPP_INCLUDED
