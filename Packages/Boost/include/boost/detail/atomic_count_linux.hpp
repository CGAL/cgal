#ifndef BOOST_DETAIL_ATOMIC_COUNT_LINUX_HPP_INCLUDED
#define BOOST_DETAIL_ATOMIC_COUNT_LINUX_HPP_INCLUDED

//
//  boost/detail/atomic_count_linux.hpp
//
//  Copyright (c) 2001, 2002 Peter Dimov and Multi Media Ltd.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

//
//  This implementation uses <asm/atomic.h>. This is a kernel header;
//  using kernel headers in a user program may cause a number of problems,
//  and not all flavors of Linux provide the atomic instructions.
//
//  This file is only provided because the performance of this implementation
//  is significantly higher than the pthreads version. Use at your own risk
//  (by defining BOOST_USE_ASM_ATOMIC_H.)
//

#include <asm/atomic.h>

namespace boost
{

namespace detail
{

class atomic_count
{
public:

    explicit atomic_count(long v)
    {
        atomic_t init = ATOMIC_INIT(v);
        value_ = init;
    }

    void operator++()
    {
        atomic_inc(&value_);
    }

    long operator--()
    {
        return !atomic_dec_and_test(&value_);
    }

    operator long() const
    {
        return atomic_read(&value_);
    }

private:

    atomic_count(atomic_count const &);
    atomic_count & operator=(atomic_count const &);

    atomic_t value_;
};

} // namespace detail

} // namespace boost

#endif // #ifndef BOOST_DETAIL_ATOMIC_COUNT_LINUX_HPP_INCLUDED
